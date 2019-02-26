% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR12
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
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR12_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR12_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:22:35
% EndTime: 2019-02-26 22:22:35
% DurationCPUTime: 0.20s
% Computational Cost: add. (356->56), mult. (536->91), div. (0->0), fcn. (676->14), ass. (0->42)
t35 = sin(qJ(3));
t38 = cos(qJ(3));
t33 = pkin(12) + qJ(5);
t30 = cos(t33);
t21 = pkin(5) * t30 + cos(pkin(12)) * pkin(4) + pkin(3);
t31 = qJ(6) + t33;
t27 = sin(t31);
t28 = cos(t31);
t42 = r_i_i_C(1) * t28 - r_i_i_C(2) * t27 + t21;
t51 = r_i_i_C(3) + pkin(11) + pkin(10) + qJ(4);
t56 = t51 * t35 + t42 * t38 + pkin(2);
t36 = sin(qJ(2));
t37 = sin(qJ(1));
t39 = cos(qJ(2));
t46 = cos(pkin(6));
t50 = cos(qJ(1));
t43 = t46 * t50;
t18 = t36 * t43 + t37 * t39;
t34 = sin(pkin(6));
t45 = t34 * t50;
t10 = t18 * t38 - t35 * t45;
t17 = t37 * t36 - t39 * t43;
t55 = (-t10 * t27 + t17 * t28) * r_i_i_C(1) + (-t10 * t28 - t17 * t27) * r_i_i_C(2);
t44 = t37 * t46;
t20 = -t36 * t44 + t50 * t39;
t49 = t34 * t37;
t14 = t20 * t38 + t35 * t49;
t19 = t50 * t36 + t39 * t44;
t5 = -t14 * t27 + t19 * t28;
t6 = t14 * t28 + t19 * t27;
t54 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t48 = t34 * t38;
t16 = t46 * t35 + t36 * t48;
t47 = t34 * t39;
t53 = (-t16 * t27 - t28 * t47) * r_i_i_C(1) + (-t16 * t28 + t27 * t47) * r_i_i_C(2);
t29 = sin(t33);
t52 = pkin(9) + pkin(5) * t29 + sin(pkin(12)) * pkin(4);
t41 = t27 * r_i_i_C(1) + t28 * r_i_i_C(2) + t52;
t9 = t18 * t35 + t38 * t45;
t15 = t34 * t36 * t35 - t46 * t38;
t13 = t20 * t35 - t37 * t48;
t1 = [-t37 * pkin(1) - t18 * pkin(2) + pkin(8) * t45 - t42 * t10 - t41 * t17 - t51 * t9, -t19 * t56 + t41 * t20, -t42 * t13 + t51 * t14, t13 (-t14 * t29 + t19 * t30) * pkin(5) + t54, t54; t50 * pkin(1) + t20 * pkin(2) + pkin(8) * t49 + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t51 * t13 + t14 * t21 + t52 * t19, -t17 * t56 + t41 * t18, t51 * t10 - t42 * t9, t9 (-t10 * t29 + t17 * t30) * pkin(5) + t55, t55; 0 (t41 * t36 + t56 * t39) * t34, -t42 * t15 + t51 * t16, t15 (-t16 * t29 - t30 * t47) * pkin(5) + t53, t53;];
Ja_transl  = t1;
