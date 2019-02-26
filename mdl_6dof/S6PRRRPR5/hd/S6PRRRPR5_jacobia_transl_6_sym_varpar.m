% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR5
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:12:54
% EndTime: 2019-02-26 20:12:55
% DurationCPUTime: 0.31s
% Computational Cost: add. (538->94), mult. (1246->168), div. (0->0), fcn. (1622->16), ass. (0->64)
t40 = sin(pkin(7));
t47 = sin(qJ(2));
t50 = cos(qJ(2));
t42 = cos(pkin(12));
t70 = cos(pkin(6));
t64 = t42 * t70;
t68 = sin(pkin(12));
t54 = t47 * t68 - t50 * t64;
t69 = cos(pkin(7));
t41 = sin(pkin(6));
t73 = t41 * t42;
t79 = t40 * t73 + t54 * t69;
t58 = t70 * t68;
t55 = t42 * t47 + t50 * t58;
t65 = t41 * t68;
t78 = -t40 * t65 + t55 * t69;
t77 = r_i_i_C(3) + pkin(11);
t76 = cos(qJ(3));
t39 = qJ(4) + pkin(13);
t37 = sin(t39);
t75 = t37 * t40;
t38 = cos(t39);
t74 = t38 * t40;
t72 = t41 * t47;
t71 = t41 * t50;
t67 = t40 * t72;
t46 = sin(qJ(3));
t63 = t46 * t69;
t62 = t70 * t40;
t45 = sin(qJ(4));
t60 = (pkin(4) * t45 + pkin(9)) * t40;
t59 = t69 * t76;
t44 = sin(qJ(6));
t48 = cos(qJ(6));
t57 = r_i_i_C(1) * t48 - r_i_i_C(2) * t44 + pkin(5);
t43 = -qJ(5) - pkin(10);
t56 = r_i_i_C(1) * t44 + r_i_i_C(2) * t48 - t43;
t49 = cos(qJ(4));
t36 = pkin(4) * t49 + pkin(3);
t53 = -t37 * t77 - t38 * t57 - t36;
t31 = t42 * t50 - t47 * t58;
t30 = t47 * t64 + t50 * t68;
t29 = -t40 * t71 + t69 * t70;
t28 = (-t47 * t63 + t50 * t76) * t41;
t27 = (t46 * t50 + t47 * t59) * t41;
t24 = t40 * t55 + t65 * t69;
t23 = t40 * t54 - t69 * t73;
t22 = t46 * t62 + (t47 * t76 + t50 * t63) * t41;
t21 = t46 * t72 - t59 * t71 - t62 * t76;
t20 = t28 * t38 + t37 * t67;
t18 = -t31 * t63 - t55 * t76;
t17 = t31 * t59 - t46 * t55;
t16 = -t30 * t63 - t54 * t76;
t15 = t30 * t59 - t46 * t54;
t14 = t31 * t76 - t46 * t78;
t13 = t31 * t46 + t76 * t78;
t12 = t30 * t76 - t46 * t79;
t11 = t30 * t46 + t76 * t79;
t10 = t22 * t38 + t29 * t37;
t8 = t18 * t38 + t31 * t75;
t6 = t16 * t38 + t30 * t75;
t4 = t14 * t38 + t24 * t37;
t2 = t12 * t38 + t23 * t37;
t1 = [0 (t17 * t44 + t48 * t8) * r_i_i_C(1) + (t17 * t48 - t44 * t8) * r_i_i_C(2) + t8 * pkin(5) + t18 * t36 - t17 * t43 - t55 * pkin(2) + t77 * (t18 * t37 - t31 * t74) + t31 * t60, t13 * t53 + t14 * t56, t77 * t4 + (-t14 * t45 + t24 * t49) * pkin(4) + t57 * (-t14 * t37 + t24 * t38) t13 (t13 * t48 - t4 * t44) * r_i_i_C(1) + (-t13 * t44 - t4 * t48) * r_i_i_C(2); 0 (t15 * t44 + t48 * t6) * r_i_i_C(1) + (t15 * t48 - t44 * t6) * r_i_i_C(2) + t6 * pkin(5) + t16 * t36 - t15 * t43 - t54 * pkin(2) + t77 * (t16 * t37 - t30 * t74) + t30 * t60, t11 * t53 + t12 * t56, t77 * t2 + (-t12 * t45 + t23 * t49) * pkin(4) + t57 * (-t12 * t37 + t23 * t38) t11 (t11 * t48 - t2 * t44) * r_i_i_C(1) + (-t11 * t44 - t2 * t48) * r_i_i_C(2); 1 (t20 * t48 + t27 * t44) * r_i_i_C(1) + (-t20 * t44 + t27 * t48) * r_i_i_C(2) + t20 * pkin(5) + t28 * t36 - t27 * t43 + t77 * (t28 * t37 - t38 * t67) + (t50 * pkin(2) + t47 * t60) * t41, t21 * t53 + t22 * t56, t77 * t10 + (-t22 * t45 + t29 * t49) * pkin(4) + t57 * (-t22 * t37 + t29 * t38) t21 (-t10 * t44 + t21 * t48) * r_i_i_C(1) + (-t10 * t48 - t21 * t44) * r_i_i_C(2);];
Ja_transl  = t1;
