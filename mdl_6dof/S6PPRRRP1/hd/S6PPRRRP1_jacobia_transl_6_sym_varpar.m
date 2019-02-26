% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRRP1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRP1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:41:33
% EndTime: 2019-02-26 19:41:33
% DurationCPUTime: 0.21s
% Computational Cost: add. (359->46), mult. (984->85), div. (0->0), fcn. (1294->14), ass. (0->50)
t68 = pkin(5) + r_i_i_C(1);
t55 = sin(pkin(12));
t56 = sin(pkin(11));
t43 = t56 * t55;
t59 = cos(pkin(12));
t60 = cos(pkin(11));
t50 = t60 * t59;
t62 = cos(pkin(6));
t35 = -t62 * t50 + t43;
t57 = sin(pkin(7));
t58 = sin(pkin(6));
t47 = t58 * t57;
t61 = cos(pkin(7));
t67 = t35 * t61 + t60 * t47;
t45 = t56 * t59;
t48 = t60 * t55;
t36 = t62 * t45 + t48;
t44 = t56 * t58;
t66 = t36 * t61 - t57 * t44;
t65 = t59 * t61 * t58 + t62 * t57;
t64 = r_i_i_C(3) + qJ(6) + pkin(10);
t63 = cos(qJ(3));
t49 = t60 * t58;
t46 = t58 * t55;
t25 = sin(qJ(5));
t28 = cos(qJ(5));
t42 = t25 * r_i_i_C(2) - t68 * t28 - pkin(4);
t41 = t28 * r_i_i_C(2) + t68 * t25 + pkin(9);
t26 = sin(qJ(4));
t29 = cos(qJ(4));
t37 = -t64 * t26 + t42 * t29 - pkin(3);
t34 = -t59 * t47 + t62 * t61;
t31 = t36 * t57 + t61 * t44;
t30 = t35 * t57 - t61 * t49;
t27 = sin(qJ(3));
t19 = -t62 * t43 + t50;
t18 = t62 * t48 + t45;
t14 = t65 * t27 + t63 * t46;
t13 = t27 * t46 - t65 * t63;
t10 = t14 * t29 + t34 * t26;
t9 = t14 * t26 - t34 * t29;
t8 = t19 * t63 - t66 * t27;
t7 = t19 * t27 + t66 * t63;
t6 = t18 * t63 - t67 * t27;
t5 = t18 * t27 + t67 * t63;
t4 = t31 * t26 + t8 * t29;
t3 = t8 * t26 - t31 * t29;
t2 = t30 * t26 + t6 * t29;
t1 = t6 * t26 - t30 * t29;
t11 = [0, t44, t37 * t7 + t41 * t8, t42 * t3 + t64 * t4 (-t7 * t25 - t4 * t28) * r_i_i_C(2) + t68 * (-t4 * t25 + t7 * t28) t3; 0, -t49, t37 * t5 + t41 * t6, t42 * t1 + t64 * t2 (-t2 * t28 - t5 * t25) * r_i_i_C(2) + t68 * (-t2 * t25 + t5 * t28) t1; 1, t62, t37 * t13 + t41 * t14, t64 * t10 + t42 * t9 (-t10 * t28 - t13 * t25) * r_i_i_C(2) + t68 * (-t10 * t25 + t13 * t28) t9;];
Ja_transl  = t11;
