% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP13
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP13_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP13_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:52:48
% EndTime: 2019-02-26 21:52:49
% DurationCPUTime: 0.19s
% Computational Cost: add. (207->51), mult. (490->80), div. (0->0), fcn. (619->10), ass. (0->37)
t45 = pkin(5) + r_i_i_C(1);
t22 = sin(qJ(5));
t26 = cos(qJ(5));
t46 = pkin(2) + pkin(9);
t47 = t26 * r_i_i_C(2) + t45 * t22 + t46;
t44 = r_i_i_C(3) + qJ(6) + pkin(10);
t20 = sin(pkin(6));
t24 = sin(qJ(2));
t43 = t20 * t24;
t25 = sin(qJ(1));
t42 = t20 * t25;
t28 = cos(qJ(2));
t41 = t20 * t28;
t29 = cos(qJ(1));
t40 = t20 * t29;
t39 = cos(pkin(6));
t38 = t20 * (pkin(3) + pkin(8));
t37 = t25 * t39;
t36 = t29 * t39;
t16 = -t24 * t37 + t29 * t28;
t15 = t29 * t24 + t28 * t37;
t23 = sin(qJ(4));
t27 = cos(qJ(4));
t4 = t15 * t23 + t27 * t42;
t1 = t16 * t26 - t4 * t22;
t19 = t26 * pkin(5) + pkin(4);
t33 = t26 * r_i_i_C(1) - t22 * r_i_i_C(2) + t19;
t13 = t25 * t24 - t28 * t36;
t7 = -t13 * t23 + t27 * t40;
t5 = t13 * t27 + t23 * t40;
t30 = t33 * t23 - t44 * t27 + qJ(3);
t14 = t24 * t36 + t25 * t28;
t12 = -t23 * t41 + t39 * t27;
t11 = t39 * t23 + t27 * t41;
t3 = -t15 * t27 + t23 * t42;
t2 = t16 * t22 + t4 * t26;
t6 = [-t25 * pkin(1) - t13 * qJ(3) - t14 * t47 + t29 * t38 + t33 * t7 + t44 * t5, -t15 * t47 + t30 * t16, t15, -t33 * t3 + t44 * t4, -t2 * r_i_i_C(2) + t1 * t45, t3; t29 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * qJ(3) + t4 * t19 + t44 * t3 + t25 * t38 + (t22 * pkin(5) + t46) * t16, -t13 * t47 + t30 * t14, t13, t33 * t5 - t44 * t7 (-t14 * t22 + t7 * t26) * r_i_i_C(2) + t45 * (t14 * t26 + t7 * t22) -t5; 0 (t30 * t24 + t47 * t28) * t20, -t41, -t33 * t11 + t44 * t12 (-t12 * t26 - t22 * t43) * r_i_i_C(2) + t45 * (-t12 * t22 + t26 * t43) t11;];
Ja_transl  = t6;
