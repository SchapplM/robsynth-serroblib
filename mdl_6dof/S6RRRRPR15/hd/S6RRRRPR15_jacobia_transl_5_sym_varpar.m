% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR15_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR15_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:57
% EndTime: 2019-02-26 22:38:57
% DurationCPUTime: 0.29s
% Computational Cost: add. (390->78), mult. (1067->135), div. (0->0), fcn. (1387->12), ass. (0->55)
t47 = sin(qJ(2));
t50 = cos(qJ(2));
t51 = cos(qJ(1));
t62 = cos(pkin(6));
t58 = t51 * t62;
t74 = sin(qJ(1));
t34 = t74 * t47 - t50 * t58;
t35 = t47 * t58 + t74 * t50;
t46 = sin(qJ(3));
t49 = cos(qJ(3));
t42 = sin(pkin(7));
t43 = sin(pkin(6));
t70 = t43 * t51;
t60 = t42 * t70;
t44 = cos(pkin(7));
t69 = t44 * t46;
t16 = t34 * t69 - t35 * t49 + t46 * t60;
t28 = -t34 * t42 + t44 * t70;
t45 = sin(qJ(4));
t48 = cos(qJ(4));
t79 = t16 * t48 + t28 * t45;
t78 = t16 * t45 - t28 * t48;
t63 = r_i_i_C(3) + qJ(5);
t77 = pkin(4) - r_i_i_C(2);
t53 = t63 * t45 + t77 * t48 + pkin(3);
t76 = pkin(11) + r_i_i_C(1);
t75 = pkin(10) * t42;
t73 = t42 * t45;
t72 = t42 * t47;
t71 = t42 * t48;
t68 = t44 * t49;
t67 = t46 * t47;
t66 = t46 * t50;
t65 = t47 * t49;
t64 = t49 * t50;
t61 = t43 * t72;
t59 = t43 * t74;
t57 = t62 * t42;
t56 = t42 * t59;
t55 = t62 * t74;
t36 = -t51 * t47 - t50 * t55;
t54 = -t36 * t42 + t44 * t59;
t52 = -t35 * t46 + (-t34 * t44 - t60) * t49;
t37 = -t47 * t55 + t51 * t50;
t33 = -t43 * t50 * t42 + t62 * t44;
t32 = (-t44 * t67 + t64) * t43;
t27 = t46 * t57 + (t44 * t66 + t65) * t43;
t22 = t36 * t49 - t37 * t69;
t20 = -t34 * t49 - t35 * t69;
t18 = t37 * t49 + (t36 * t44 + t56) * t46;
t17 = -t36 * t68 + t37 * t46 - t49 * t56;
t11 = t27 * t45 - t33 * t48;
t6 = t18 * t48 + t54 * t45;
t5 = t18 * t45 - t54 * t48;
t1 = [-t74 * pkin(1) - t35 * pkin(2) + t16 * pkin(3) + pkin(9) * t70 + t28 * pkin(10) + t76 * t52 + t63 * t78 + t77 * t79, t37 * t75 + t36 * pkin(2) + t22 * pkin(3) + t63 * (t22 * t45 - t37 * t71) + t76 * (t36 * t46 + t37 * t68) + t77 * (t22 * t48 + t37 * t73) -t53 * t17 + t76 * t18, -t77 * t5 + t63 * t6, t5, 0; t51 * pkin(1) + t37 * pkin(2) + t18 * pkin(3) + pkin(9) * t59 + t54 * pkin(10) + t76 * t17 + t63 * t5 + t77 * t6, t35 * t75 - t34 * pkin(2) + t20 * pkin(3) + t77 * (t20 * t48 + t35 * t73) + t63 * (t20 * t45 - t35 * t71) + t76 * (-t34 * t46 + t35 * t68) -t16 * t76 + t53 * t52, -t63 * t79 + t77 * t78, -t78, 0; 0, t32 * pkin(3) + t77 * (t32 * t48 + t45 * t61) + t63 * (t32 * t45 - t48 * t61) + (pkin(2) * t50 + pkin(10) * t72 + t76 * (t44 * t65 + t66)) * t43, t76 * t27 + t53 * (t49 * t57 + (t44 * t64 - t67) * t43) t63 * (t27 * t48 + t33 * t45) - t77 * t11, t11, 0;];
Ja_transl  = t1;
