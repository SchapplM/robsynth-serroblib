% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:25
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR14_jacobia_transl_4_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_transl_4_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobia_transl_4_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_transl_4_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:25:33
% EndTime: 2019-01-03 10:25:33
% DurationCPUTime: 0.39s
% Computational Cost: add. (273->80), mult. (775->145), div. (0->0), fcn. (1007->14), ass. (0->56)
t34 = sin(qJ(4));
t37 = cos(qJ(4));
t35 = sin(qJ(2));
t36 = sin(qJ(1));
t38 = cos(qJ(2));
t39 = cos(qJ(1));
t56 = cos(pkin(6));
t54 = t39 * t56;
t21 = t36 * t35 - t38 * t54;
t29 = sin(pkin(7));
t33 = cos(pkin(7));
t30 = sin(pkin(6));
t61 = t30 * t39;
t15 = -t21 * t29 + t33 * t61;
t28 = sin(pkin(8));
t32 = cos(pkin(8));
t22 = t35 * t54 + t36 * t38;
t27 = sin(pkin(14));
t31 = cos(pkin(14));
t41 = t21 * t33 + t29 * t61;
t5 = t22 * t27 + t31 * t41;
t50 = t15 * t28 + t32 * t5;
t6 = -t22 * t31 + t27 * t41;
t73 = t50 * t34 + t6 * t37;
t72 = -t6 * t34 + t50 * t37;
t68 = pkin(11) + r_i_i_C(3);
t65 = t27 * t33;
t64 = t28 * t29;
t63 = t29 * t32;
t62 = t30 * t36;
t60 = t31 * t33;
t59 = t33 * t35;
t58 = t33 * t38;
t57 = t29 * qJ(3);
t55 = t36 * t56;
t53 = t56 * t29;
t52 = t34 * r_i_i_C(1) + t37 * r_i_i_C(2);
t23 = -t39 * t35 - t38 * t55;
t17 = -t23 * t29 + t33 * t62;
t24 = -t35 * t55 + t39 * t38;
t40 = t23 * t33 + t29 * t62;
t7 = -t24 * t27 + t31 * t40;
t47 = t17 * t28 + t32 * t7;
t20 = -t30 * t38 * t29 + t56 * t33;
t46 = (t31 * t53 + (-t27 * t35 + t31 * t58) * t30) * t32 + t20 * t28;
t9 = t21 * t27 - t22 * t60;
t44 = t22 * t64 + t32 * t9;
t11 = -t23 * t27 - t24 * t60;
t42 = t11 * t32 + t24 * t64;
t14 = t30 * t35 * t31 + (t30 * t58 + t53) * t27;
t12 = t23 * t31 - t24 * t65;
t10 = -t21 * t31 - t22 * t65;
t8 = t24 * t31 + t27 * t40;
t2 = t47 * t34 + t8 * t37;
t1 = -t8 * t34 + t47 * t37;
t3 = [t73 * r_i_i_C(1) + t72 * r_i_i_C(2) + t6 * pkin(3) - t22 * pkin(2) - t36 * pkin(1) + pkin(10) * t61 + t15 * qJ(3) + t68 * (t15 * t32 - t5 * t28) (t12 * t37 + t42 * t34) * r_i_i_C(1) + (-t12 * t34 + t42 * t37) * r_i_i_C(2) + t12 * pkin(3) + t23 * pkin(2) + t24 * t57 + t68 * (-t11 * t28 + t24 * t63) t17, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t39 * pkin(1) + t24 * pkin(2) + t8 * pkin(3) + pkin(10) * t62 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t17 * qJ(3) + t68 * (t17 * t32 - t7 * t28) (t10 * t37 + t44 * t34) * r_i_i_C(1) + (-t10 * t34 + t44 * t37) * r_i_i_C(2) + t10 * pkin(3) - t21 * pkin(2) + t22 * t57 + t68 * (t22 * t63 - t9 * t28) -t15, -t72 * r_i_i_C(1) + t73 * r_i_i_C(2), 0, 0; 0 ((t37 * r_i_i_C(1) - t34 * r_i_i_C(2) + pkin(3)) * (-t27 * t59 + t31 * t38) + (-t68 * t28 + t52 * t32) * (-t27 * t38 - t31 * t59) + t38 * pkin(2) + (t52 * t28 + t68 * t32 + qJ(3)) * t35 * t29) * t30, t20 (-t14 * t34 + t46 * t37) * r_i_i_C(1) + (-t14 * t37 - t46 * t34) * r_i_i_C(2), 0, 0;];
Ja_transl  = t3;
