% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:46:00
% EndTime: 2019-02-26 21:46:00
% DurationCPUTime: 0.12s
% Computational Cost: add. (195->39), mult. (152->46), div. (0->0), fcn. (165->10), ass. (0->32)
t21 = qJ(2) + pkin(10);
t18 = qJ(4) + t21;
t15 = sin(t18);
t16 = cos(t18);
t23 = sin(qJ(5));
t42 = r_i_i_C(2) * t23;
t47 = r_i_i_C(3) * t16 + t15 * t42;
t25 = cos(qJ(5));
t17 = pkin(5) * t25 + pkin(4);
t22 = -qJ(6) - pkin(9);
t46 = t16 * t17 + (r_i_i_C(3) - t22) * t15;
t45 = pkin(5) + r_i_i_C(1);
t33 = pkin(3) * cos(t21) + cos(qJ(2)) * pkin(2);
t44 = pkin(1) + t33 + t46;
t43 = r_i_i_C(1) * t25;
t24 = sin(qJ(1));
t39 = t47 * t24;
t26 = cos(qJ(1));
t38 = t47 * t26;
t37 = t23 * t26;
t36 = t24 * t23;
t35 = t24 * t25;
t34 = t25 * t26;
t31 = pkin(5) * t23 + pkin(7) + pkin(8) + qJ(3);
t3 = -t16 * t37 + t35;
t1 = t16 * t36 + t34;
t29 = -t16 * t22 + (-t17 - t43) * t15;
t28 = (-t42 + t43) * t16 + t46;
t27 = -pkin(3) * sin(t21) - sin(qJ(2)) * pkin(2) + t29;
t4 = t16 * t34 + t36;
t2 = -t16 * t35 + t37;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t44 * t24 + t31 * t26, t27 * t26 + t38, t24, t29 * t26 + t38, -t4 * r_i_i_C(2) + t45 * t3, t26 * t15; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t31 * t24 + t44 * t26, t27 * t24 + t39, -t26, t29 * t24 + t39, t2 * r_i_i_C(2) - t45 * t1, t24 * t15; 0, t28 + t33, 0, t28 (-r_i_i_C(2) * t25 - t45 * t23) * t15, -t16;];
Ja_transl  = t5;
