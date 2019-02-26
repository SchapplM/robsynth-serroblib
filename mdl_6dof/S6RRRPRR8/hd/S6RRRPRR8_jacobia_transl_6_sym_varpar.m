% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:57
% EndTime: 2019-02-26 22:19:57
% DurationCPUTime: 0.22s
% Computational Cost: add. (379->64), mult. (547->100), div. (0->0), fcn. (686->14), ass. (0->45)
t39 = cos(qJ(3));
t25 = t39 * pkin(3) + pkin(2);
t30 = qJ(3) + pkin(12);
t26 = sin(t30);
t27 = cos(t30);
t38 = cos(qJ(5));
t24 = t38 * pkin(5) + pkin(4);
t31 = qJ(5) + qJ(6);
t28 = sin(t31);
t29 = cos(t31);
t46 = r_i_i_C(1) * t29 - r_i_i_C(2) * t28 + t24;
t56 = r_i_i_C(3) + pkin(11) + pkin(10);
t60 = t56 * t26 + t46 * t27 + t25;
t36 = sin(qJ(2));
t37 = sin(qJ(1));
t40 = cos(qJ(2));
t41 = cos(qJ(1));
t51 = cos(pkin(6));
t48 = t41 * t51;
t18 = t36 * t48 + t37 * t40;
t32 = sin(pkin(6));
t52 = t32 * t41;
t10 = t18 * t27 - t26 * t52;
t17 = t37 * t36 - t40 * t48;
t59 = (-t10 * t28 + t17 * t29) * r_i_i_C(1) + (-t10 * t29 - t17 * t28) * r_i_i_C(2);
t49 = t37 * t51;
t20 = -t36 * t49 + t41 * t40;
t54 = t32 * t37;
t14 = t20 * t27 + t26 * t54;
t19 = t41 * t36 + t40 * t49;
t5 = -t14 * t28 + t19 * t29;
t6 = t14 * t29 + t19 * t28;
t58 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t55 = t32 * t36;
t16 = t51 * t26 + t27 * t55;
t53 = t32 * t40;
t57 = (-t16 * t28 - t29 * t53) * r_i_i_C(1) + (-t16 * t29 + t28 * t53) * r_i_i_C(2);
t34 = sin(qJ(5));
t50 = t34 * pkin(5) + pkin(9) + qJ(4);
t35 = sin(qJ(3));
t47 = t32 * (pkin(3) * t35 + pkin(8));
t45 = -t18 * t26 - t27 * t52;
t44 = t28 * r_i_i_C(1) + t29 * r_i_i_C(2) + t50;
t13 = t20 * t26 - t27 * t54;
t1 = [-t37 * pkin(1) - t46 * t10 - t44 * t17 - t18 * t25 + t41 * t47 + t56 * t45, -t19 * t60 + t44 * t20, t56 * t14 + (-t20 * t35 + t39 * t54) * pkin(3) - t46 * t13, t19 (-t14 * t34 + t19 * t38) * pkin(5) + t58, t58; t41 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t56 * t13 + t14 * t24 + t50 * t19 + t20 * t25 + t37 * t47, -t17 * t60 + t44 * t18, t56 * t10 + (-t18 * t35 - t39 * t52) * pkin(3) + t46 * t45, t17 (-t10 * t34 + t17 * t38) * pkin(5) + t59, t59; 0 (t44 * t36 + t60 * t40) * t32, t56 * t16 + (-t35 * t55 + t51 * t39) * pkin(3) + t46 * (-t26 * t55 + t51 * t27) -t53 (-t16 * t34 - t38 * t53) * pkin(5) + t57, t57;];
Ja_transl  = t1;
