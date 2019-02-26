% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRRR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobia_transl_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:55
% EndTime: 2019-02-26 19:43:56
% DurationCPUTime: 0.41s
% Computational Cost: add. (1073->96), mult. (3094->175), div. (0->0), fcn. (4112->18), ass. (0->74)
t90 = r_i_i_C(3) + pkin(12);
t42 = sin(pkin(8));
t41 = sin(pkin(13));
t43 = cos(pkin(14));
t44 = cos(pkin(6));
t40 = sin(pkin(14));
t82 = cos(pkin(13));
t77 = t82 * t40;
t37 = t41 * t43 + t44 * t77;
t48 = sin(qJ(3));
t51 = cos(qJ(3));
t76 = t82 * t43;
t70 = -t41 * t40 + t44 * t76;
t80 = sin(pkin(7));
t81 = sin(pkin(6));
t73 = t80 * t81;
t84 = cos(pkin(7));
t91 = t70 * t84 - t82 * t73;
t54 = t37 * t48 - t91 * t51;
t74 = t81 * t84;
t58 = -t70 * t80 - t82 * t74;
t83 = cos(pkin(8));
t94 = -t58 * t42 + t54 * t83;
t87 = t41 * t44;
t38 = -t40 * t87 + t76;
t69 = -t43 * t87 - t77;
t62 = t41 * t73 + t69 * t84;
t55 = t38 * t48 - t62 * t51;
t61 = t41 * t74 - t69 * t80;
t93 = -t61 * t42 + t55 * t83;
t65 = t43 * t74 + t80 * t44;
t79 = t40 * t81;
t60 = t48 * t79 - t65 * t51;
t66 = -t43 * t73 + t44 * t84;
t92 = -t66 * t42 + t60 * t83;
t89 = pkin(10) * t42;
t88 = cos(qJ(4));
t46 = sin(qJ(5));
t86 = t42 * t46;
t50 = cos(qJ(5));
t85 = t42 * t50;
t47 = sin(qJ(4));
t78 = t47 * t83;
t75 = t83 * t88;
t45 = sin(qJ(6));
t49 = cos(qJ(6));
t72 = t49 * r_i_i_C(1) - t45 * r_i_i_C(2) + pkin(5);
t71 = t45 * r_i_i_C(1) + t49 * r_i_i_C(2) + pkin(11);
t64 = -t90 * t46 - t72 * t50 - pkin(4);
t35 = t65 * t48 + t51 * t79;
t31 = t60 * t42 + t66 * t83;
t30 = t38 * t51 + t62 * t48;
t29 = t37 * t51 + t91 * t48;
t26 = -t35 * t78 - t60 * t88;
t25 = t35 * t75 - t60 * t47;
t24 = t55 * t42 + t61 * t83;
t23 = t54 * t42 + t58 * t83;
t22 = t35 * t88 - t92 * t47;
t21 = t35 * t47 + t92 * t88;
t20 = t26 * t50 + t35 * t86;
t18 = -t30 * t78 - t55 * t88;
t17 = t30 * t75 - t55 * t47;
t16 = -t29 * t78 - t54 * t88;
t15 = t29 * t75 - t54 * t47;
t14 = t22 * t50 + t31 * t46;
t12 = t30 * t88 - t93 * t47;
t11 = t30 * t47 + t93 * t88;
t10 = t29 * t88 - t94 * t47;
t9 = t29 * t47 + t94 * t88;
t8 = t18 * t50 + t30 * t86;
t6 = t16 * t50 + t29 * t86;
t4 = t12 * t50 + t24 * t46;
t2 = t10 * t50 + t23 * t46;
t1 = [0, t41 * t81 (t17 * t45 + t8 * t49) * r_i_i_C(1) + (t17 * t49 - t8 * t45) * r_i_i_C(2) + t8 * pkin(5) + t18 * pkin(4) + t17 * pkin(11) - t55 * pkin(3) + t30 * t89 + t90 * (t18 * t46 - t30 * t85) t64 * t11 + t71 * t12, t90 * t4 + t72 * (-t12 * t46 + t24 * t50) (t11 * t49 - t4 * t45) * r_i_i_C(1) + (-t11 * t45 - t4 * t49) * r_i_i_C(2); 0, -t82 * t81 (t15 * t45 + t6 * t49) * r_i_i_C(1) + (t15 * t49 - t6 * t45) * r_i_i_C(2) + t6 * pkin(5) + t16 * pkin(4) + t15 * pkin(11) - t54 * pkin(3) + t29 * t89 + t90 * (t16 * t46 - t29 * t85) t71 * t10 + t64 * t9, t90 * t2 + t72 * (-t10 * t46 + t23 * t50) (-t2 * t45 + t9 * t49) * r_i_i_C(1) + (-t2 * t49 - t9 * t45) * r_i_i_C(2); 1, t44 (t20 * t49 + t25 * t45) * r_i_i_C(1) + (-t20 * t45 + t25 * t49) * r_i_i_C(2) + t20 * pkin(5) + t26 * pkin(4) + t25 * pkin(11) - t60 * pkin(3) + t35 * t89 + t90 * (t26 * t46 - t35 * t85) t64 * t21 + t71 * t22, t90 * t14 + t72 * (-t22 * t46 + t31 * t50) (-t14 * t45 + t21 * t49) * r_i_i_C(1) + (-t14 * t49 - t21 * t45) * r_i_i_C(2);];
Ja_transl  = t1;
