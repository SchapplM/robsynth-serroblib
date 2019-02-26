% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR11_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR11_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:49
% EndTime: 2019-02-26 22:21:50
% DurationCPUTime: 0.28s
% Computational Cost: add. (409->67), mult. (1057->107), div. (0->0), fcn. (1379->12), ass. (0->47)
t43 = sin(qJ(5));
t44 = sin(qJ(3));
t48 = cos(qJ(5));
t49 = cos(qJ(3));
t42 = sin(qJ(6));
t47 = cos(qJ(6));
t55 = t47 * r_i_i_C(1) - t42 * r_i_i_C(2) + pkin(5);
t72 = r_i_i_C(3) + pkin(11);
t81 = t72 * (t43 * t49 - t44 * t48) + (t43 * t44 + t48 * t49) * t55;
t45 = sin(qJ(2));
t46 = sin(qJ(1));
t50 = cos(qJ(2));
t67 = cos(pkin(6));
t71 = cos(qJ(1));
t60 = t67 * t71;
t33 = t45 * t60 + t46 * t50;
t41 = sin(pkin(6));
t63 = t41 * t71;
t22 = t33 * t44 + t49 * t63;
t23 = t33 * t49 - t44 * t63;
t80 = t22 * t48 - t23 * t43;
t6 = t22 * t43 + t23 * t48;
t74 = pkin(3) + pkin(4);
t77 = t44 * qJ(4) + t74 * t49 + pkin(2);
t79 = -t77 - t81;
t61 = t46 * t67;
t35 = -t45 * t61 + t71 * t50;
t69 = t41 * t49;
t26 = t35 * t44 - t46 * t69;
t70 = t41 * t46;
t27 = t35 * t49 + t44 * t70;
t11 = -t26 * t48 + t27 * t43;
t12 = t26 * t43 + t27 * t48;
t78 = -t55 * t11 + t72 * t12;
t30 = t41 * t45 * t44 - t67 * t49;
t31 = t67 * t44 + t45 * t69;
t20 = t30 * t43 + t31 * t48;
t76 = t55 * (t30 * t48 - t31 * t43) + t72 * t20;
t75 = t55 * t80 + t72 * t6;
t73 = pkin(9) - pkin(10);
t68 = t41 * t50;
t54 = t42 * r_i_i_C(1) + t47 * r_i_i_C(2) - t73;
t34 = t71 * t45 + t50 * t61;
t32 = t46 * t45 - t50 * t60;
t2 = t12 * t47 - t34 * t42;
t1 = -t12 * t42 - t34 * t47;
t3 = [-t46 * pkin(1) - t33 * pkin(2) + pkin(8) * t63 - t22 * qJ(4) - t74 * t23 + t54 * t32 - t55 * t6 + t72 * t80, t79 * t34 - t35 * t54, t27 * qJ(4) - t74 * t26 - t78, t26, t78, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t71 * pkin(1) + t35 * pkin(2) + t12 * pkin(5) + pkin(8) * t70 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t26 * qJ(4) + t72 * t11 + t74 * t27 + t73 * t34, t79 * t32 - t33 * t54, t23 * qJ(4) - t74 * t22 - t75, t22, t75 (-t32 * t47 - t6 * t42) * r_i_i_C(1) + (t32 * t42 - t6 * t47) * r_i_i_C(2); 0 (-t54 * t45 + t77 * t50) * t41 + t81 * t68, t31 * qJ(4) - t74 * t30 - t76, t30, t76 (-t20 * t42 + t47 * t68) * r_i_i_C(1) + (-t20 * t47 - t42 * t68) * r_i_i_C(2);];
Ja_transl  = t3;
