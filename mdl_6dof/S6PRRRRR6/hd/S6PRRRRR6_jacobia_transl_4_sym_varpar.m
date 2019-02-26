% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRRR6
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRR6_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR6_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobia_transl_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:51
% EndTime: 2019-02-26 20:21:51
% DurationCPUTime: 0.30s
% Computational Cost: add. (250->71), mult. (730->136), div. (0->0), fcn. (944->14), ass. (0->57)
t26 = cos(pkin(8));
t28 = sin(qJ(4));
t31 = cos(qJ(4));
t45 = r_i_i_C(1) * t28 + r_i_i_C(2) * t31;
t22 = sin(pkin(8));
t65 = pkin(11) + r_i_i_C(3);
t49 = t65 * t22;
t34 = -t26 * t45 + t49;
t23 = sin(pkin(7));
t64 = t23 * pkin(10);
t25 = cos(pkin(14));
t30 = sin(qJ(2));
t33 = cos(qJ(2));
t51 = sin(pkin(14));
t52 = cos(pkin(6));
t41 = t52 * t51;
t19 = t25 * t33 - t30 * t41;
t29 = sin(qJ(3));
t32 = cos(qJ(3));
t18 = -t25 * t30 - t33 * t41;
t27 = cos(pkin(7));
t24 = sin(pkin(6));
t48 = t24 * t51;
t35 = t18 * t27 + t23 * t48;
t4 = t19 * t32 + t29 * t35;
t63 = t4 * t28;
t62 = t4 * t31;
t61 = t22 * t23;
t60 = t23 * t26;
t59 = t24 * t25;
t58 = t27 * t29;
t57 = t27 * t32;
t56 = t29 * t30;
t55 = t29 * t33;
t54 = t30 * t32;
t53 = t32 * t33;
t50 = t23 * t59;
t47 = t25 * t52;
t46 = t52 * t23;
t16 = -t30 * t51 + t33 * t47;
t17 = t30 * t47 + t33 * t51;
t1 = -t17 * t29 + (t16 * t27 - t50) * t32;
t44 = t1 * t26 + (-t16 * t23 - t27 * t59) * t22;
t3 = -t19 * t29 + t32 * t35;
t43 = (-t18 * t23 + t27 * t48) * t22 + t26 * t3;
t9 = t32 * t46 + (t27 * t53 - t56) * t24;
t42 = (-t23 * t24 * t33 + t27 * t52) * t22 + t26 * t9;
t40 = r_i_i_C(1) * t31 - r_i_i_C(2) * t28 + pkin(3);
t5 = -t16 * t29 - t17 * t57;
t38 = t17 * t61 + t26 * t5;
t7 = -t18 * t29 - t19 * t57;
t36 = t19 * t61 + t26 * t7;
t10 = t29 * t46 + (t27 * t55 + t54) * t24;
t8 = t18 * t32 - t19 * t58;
t6 = t16 * t32 - t17 * t58;
t2 = t16 * t58 + t17 * t32 - t29 * t50;
t11 = [0 (t28 * t36 + t8 * t31) * r_i_i_C(1) + (-t8 * t28 + t31 * t36) * r_i_i_C(2) + t8 * pkin(3) + t18 * pkin(2) + t19 * t64 + t65 * (t19 * t60 - t22 * t7) (-t26 * t63 + t3 * t31) * r_i_i_C(1) + (-t26 * t62 - t28 * t3) * r_i_i_C(2) + t3 * pkin(3) + t4 * t49 (t31 * t43 - t63) * r_i_i_C(1) + (-t28 * t43 - t62) * r_i_i_C(2), 0, 0; 0 (t28 * t38 + t6 * t31) * r_i_i_C(1) + (-t6 * t28 + t31 * t38) * r_i_i_C(2) + t6 * pkin(3) + t16 * pkin(2) + t17 * t64 + t65 * (t17 * t60 - t22 * t5) t1 * t40 + t2 * t34 (-t2 * t28 + t31 * t44) * r_i_i_C(1) + (-t2 * t31 - t28 * t44) * r_i_i_C(2), 0, 0; 1 (t40 * (-t27 * t56 + t53) - t34 * (-t27 * t54 - t55) + t33 * pkin(2) + (t22 * t45 + t26 * t65 + pkin(10)) * t30 * t23) * t24, t10 * t34 + t40 * t9 (-t10 * t28 + t31 * t42) * r_i_i_C(1) + (-t10 * t31 - t28 * t42) * r_i_i_C(2), 0, 0;];
Ja_transl  = t11;
