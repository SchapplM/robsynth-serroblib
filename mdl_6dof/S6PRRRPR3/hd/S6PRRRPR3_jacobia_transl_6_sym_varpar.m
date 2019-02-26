% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:42
% EndTime: 2019-02-26 20:11:43
% DurationCPUTime: 0.19s
% Computational Cost: add. (292->42), mult. (463->73), div. (0->0), fcn. (581->12), ass. (0->34)
t59 = pkin(4) + pkin(10) + r_i_i_C(3);
t33 = sin(qJ(6));
t36 = cos(qJ(6));
t58 = t33 * r_i_i_C(1) + t36 * r_i_i_C(2) + qJ(5);
t30 = qJ(3) + qJ(4);
t28 = sin(t30);
t29 = cos(t30);
t37 = cos(qJ(3));
t57 = t37 * pkin(3) + t58 * t28 + t59 * t29 + pkin(2);
t31 = sin(pkin(11));
t32 = sin(pkin(6));
t54 = t31 * t32;
t35 = sin(qJ(2));
t53 = t32 * t35;
t38 = cos(qJ(2));
t52 = t32 * t38;
t51 = cos(pkin(6));
t50 = cos(pkin(11));
t48 = t31 * t51;
t47 = t32 * t50;
t46 = t51 * t50;
t44 = t36 * r_i_i_C(1) - t33 * r_i_i_C(2) + pkin(5) + pkin(8) + pkin(9);
t20 = t31 * t38 + t35 * t46;
t9 = t20 * t28 + t29 * t47;
t43 = -t59 * t9 + t58 * (t20 * t29 - t28 * t47);
t22 = -t35 * t48 + t50 * t38;
t11 = t22 * t28 - t29 * t54;
t42 = t58 * (t22 * t29 + t28 * t54) - t59 * t11;
t17 = t28 * t53 - t51 * t29;
t41 = t58 * (t51 * t28 + t29 * t53) - t59 * t17;
t34 = sin(qJ(3));
t21 = t50 * t35 + t38 * t48;
t19 = t31 * t35 - t38 * t46;
t1 = [0, -t21 * t57 + t44 * t22 (-t22 * t34 + t37 * t54) * pkin(3) + t42, t42, t11 (t11 * t36 - t21 * t33) * r_i_i_C(1) + (-t11 * t33 - t21 * t36) * r_i_i_C(2); 0, -t19 * t57 + t44 * t20 (-t20 * t34 - t37 * t47) * pkin(3) + t43, t43, t9 (-t19 * t33 + t9 * t36) * r_i_i_C(1) + (-t19 * t36 - t9 * t33) * r_i_i_C(2); 1 (t44 * t35 + t57 * t38) * t32 (-t34 * t53 + t51 * t37) * pkin(3) + t41, t41, t17 (t17 * t36 + t33 * t52) * r_i_i_C(1) + (-t17 * t33 + t36 * t52) * r_i_i_C(2);];
Ja_transl  = t1;
