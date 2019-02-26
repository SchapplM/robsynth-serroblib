% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRP5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:17:28
% EndTime: 2019-02-26 20:17:28
% DurationCPUTime: 0.32s
% Computational Cost: add. (496->73), mult. (1357->138), div. (0->0), fcn. (1769->14), ass. (0->57)
t83 = pkin(5) + r_i_i_C(1);
t37 = sin(pkin(7));
t44 = sin(qJ(2));
t47 = cos(qJ(2));
t39 = cos(pkin(12));
t72 = cos(pkin(6));
t66 = t39 * t72;
t70 = sin(pkin(12));
t54 = t44 * t70 - t47 * t66;
t71 = cos(pkin(7));
t38 = sin(pkin(6));
t75 = t38 * t39;
t82 = t37 * t75 + t54 * t71;
t58 = t72 * t70;
t55 = t39 * t44 + t47 * t58;
t67 = t38 * t70;
t81 = -t37 * t67 + t55 * t71;
t80 = pkin(9) * t37;
t79 = r_i_i_C(3) + qJ(6) + pkin(11);
t78 = cos(qJ(3));
t42 = sin(qJ(4));
t77 = t37 * t42;
t46 = cos(qJ(4));
t76 = t37 * t46;
t74 = t38 * t44;
t73 = t38 * t47;
t69 = t37 * t74;
t43 = sin(qJ(3));
t65 = t43 * t71;
t64 = t72 * t37;
t60 = t71 * t78;
t41 = sin(qJ(5));
t45 = cos(qJ(5));
t57 = -t41 * r_i_i_C(2) + t45 * t83 + pkin(4);
t56 = t45 * r_i_i_C(2) + t41 * t83 + pkin(10);
t53 = -t37 * t73 + t71 * t72;
t52 = -t42 * t79 - t46 * t57 - pkin(3);
t49 = t37 * t54 - t71 * t75;
t48 = t37 * t55 + t67 * t71;
t31 = t39 * t47 - t44 * t58;
t30 = t44 * t66 + t47 * t70;
t28 = (-t44 * t65 + t47 * t78) * t38;
t24 = t43 * t64 + (t44 * t78 + t47 * t65) * t38;
t23 = t43 * t74 - t60 * t73 - t64 * t78;
t18 = -t31 * t65 - t55 * t78;
t16 = -t30 * t65 - t54 * t78;
t14 = t24 * t46 + t42 * t53;
t13 = t24 * t42 - t46 * t53;
t12 = t31 * t78 - t43 * t81;
t11 = t31 * t43 + t78 * t81;
t10 = t30 * t78 - t43 * t82;
t9 = t30 * t43 + t78 * t82;
t4 = t12 * t46 + t42 * t48;
t3 = t12 * t42 - t46 * t48;
t2 = t10 * t46 + t42 * t49;
t1 = t10 * t42 - t46 * t49;
t5 = [0, t18 * pkin(3) - t55 * pkin(2) + t31 * t80 + t79 * (t18 * t42 - t31 * t76) + t57 * (t18 * t46 + t31 * t77) + t56 * (t31 * t60 - t43 * t55) t11 * t52 + t12 * t56, -t3 * t57 + t4 * t79 (-t11 * t41 - t4 * t45) * r_i_i_C(2) + t83 * (t11 * t45 - t4 * t41) t3; 0, t16 * pkin(3) - t54 * pkin(2) + t30 * t80 + t79 * (t16 * t42 - t30 * t76) + t57 * (t16 * t46 + t30 * t77) + t56 * (t30 * t60 - t43 * t54) t10 * t56 + t52 * t9, -t1 * t57 + t2 * t79 (-t2 * t45 - t41 * t9) * r_i_i_C(2) + t83 * (-t2 * t41 + t45 * t9) t1; 1, t28 * pkin(3) + t57 * (t28 * t46 + t42 * t69) + t79 * (t28 * t42 - t46 * t69) + (t47 * pkin(2) + t44 * t80 + t56 * (t43 * t47 + t44 * t60)) * t38, t23 * t52 + t24 * t56, -t13 * t57 + t14 * t79 (-t14 * t45 - t23 * t41) * r_i_i_C(2) + t83 * (-t14 * t41 + t23 * t45) t13;];
Ja_transl  = t5;
