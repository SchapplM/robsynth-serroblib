% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:07
% EndTime: 2019-02-26 20:08:07
% DurationCPUTime: 0.26s
% Computational Cost: add. (441->76), mult. (1239->140), div. (0->0), fcn. (1617->14), ass. (0->56)
t73 = r_i_i_C(3) + pkin(11);
t38 = sin(pkin(7));
t39 = sin(pkin(6));
t72 = t38 * t39;
t43 = sin(qJ(5));
t71 = t38 * t43;
t47 = cos(qJ(5));
t70 = t38 * t47;
t41 = cos(pkin(7));
t69 = t39 * t41;
t44 = sin(qJ(3));
t68 = t41 * t44;
t48 = cos(qJ(3));
t67 = t41 * t48;
t45 = sin(qJ(2));
t66 = t44 * t45;
t49 = cos(qJ(2));
t65 = t44 * t49;
t64 = t45 * t48;
t63 = t48 * t49;
t62 = cos(pkin(6));
t61 = sin(pkin(12));
t40 = cos(pkin(12));
t60 = t40 * t72;
t59 = t45 * t72;
t58 = (pkin(4) + pkin(9)) * t38;
t57 = t39 * t61;
t56 = t40 * t62;
t55 = t62 * t38;
t54 = t38 * t57;
t53 = t62 * t61;
t42 = sin(qJ(6));
t46 = cos(qJ(6));
t52 = t46 * r_i_i_C(1) - t42 * r_i_i_C(2) + pkin(5);
t51 = t42 * r_i_i_C(1) + t46 * r_i_i_C(2) + pkin(3) + pkin(10);
t50 = t52 * t43 - t73 * t47 + qJ(4);
t33 = t40 * t49 - t45 * t53;
t32 = -t40 * t45 - t49 * t53;
t31 = t45 * t56 + t61 * t49;
t30 = -t61 * t45 + t49 * t56;
t29 = t62 * t41 - t49 * t72;
t27 = (t41 * t64 + t65) * t39;
t24 = -t32 * t38 + t41 * t57;
t23 = -t30 * t38 - t40 * t69;
t22 = t44 * t55 + (t41 * t65 + t64) * t39;
t21 = t39 * t66 - t48 * t55 - t63 * t69;
t17 = t32 * t44 + t33 * t67;
t15 = t30 * t44 + t31 * t67;
t14 = t21 * t43 + t29 * t47;
t12 = t33 * t48 + (t32 * t41 + t54) * t44;
t11 = -t32 * t67 + t33 * t44 - t48 * t54;
t10 = t31 * t48 + (t30 * t41 - t60) * t44;
t9 = -t30 * t67 + t31 * t44 + t48 * t60;
t4 = t11 * t43 + t24 * t47;
t2 = t23 * t47 + t9 * t43;
t1 = [0, t32 * pkin(2) + t17 * qJ(4) - t73 * (t17 * t47 - t33 * t71) + t33 * t58 + t52 * (t17 * t43 + t33 * t70) + t51 * (t32 * t48 - t33 * t68) -t51 * t11 + t50 * t12, t11, t73 * t4 + t52 * (t11 * t47 - t24 * t43) (t12 * t46 - t4 * t42) * r_i_i_C(1) + (-t12 * t42 - t4 * t46) * r_i_i_C(2); 0, t30 * pkin(2) + t15 * qJ(4) - t73 * (t15 * t47 - t31 * t71) + t31 * t58 + t52 * (t15 * t43 + t31 * t70) + t51 * (t30 * t48 - t31 * t68) t50 * t10 - t51 * t9, t9, t73 * t2 + t52 * (-t23 * t43 + t9 * t47) (t10 * t46 - t2 * t42) * r_i_i_C(1) + (-t10 * t42 - t2 * t46) * r_i_i_C(2); 1, t27 * qJ(4) + t52 * (t27 * t43 + t47 * t59) + t73 * (-t27 * t47 + t43 * t59) + (t51 * (-t41 * t66 + t63) + t49 * pkin(2) + t45 * t58) * t39, -t51 * t21 + t50 * t22, t21, t73 * t14 + t52 * (t21 * t47 - t29 * t43) (-t14 * t42 + t22 * t46) * r_i_i_C(1) + (-t14 * t46 - t22 * t42) * r_i_i_C(2);];
Ja_transl  = t1;
