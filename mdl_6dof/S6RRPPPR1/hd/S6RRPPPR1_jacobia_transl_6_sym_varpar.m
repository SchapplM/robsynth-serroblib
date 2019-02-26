% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPPR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:48
% EndTime: 2019-02-26 21:21:49
% DurationCPUTime: 0.15s
% Computational Cost: add. (169->39), mult. (224->55), div. (0->0), fcn. (272->10), ass. (0->29)
t26 = -r_i_i_C(3) - pkin(8) + qJ(4);
t12 = qJ(2) + pkin(9);
t9 = sin(t12);
t34 = cos(qJ(2)) * pkin(2) + t26 * t9;
t10 = cos(t12);
t33 = t10 * pkin(3) + pkin(1) + t34;
t13 = sin(pkin(10));
t14 = cos(pkin(10));
t16 = sin(qJ(6));
t19 = cos(qJ(6));
t31 = pkin(4) + pkin(5);
t22 = t19 * r_i_i_C(1) - t16 * r_i_i_C(2) + t31;
t23 = t16 * r_i_i_C(1) + t19 * r_i_i_C(2) + qJ(5);
t32 = t23 * t13 + t22 * t14 + pkin(3);
t18 = sin(qJ(1));
t30 = t18 * t13;
t29 = t18 * t14;
t20 = cos(qJ(1));
t28 = t20 * t13;
t27 = t20 * t14;
t21 = -sin(qJ(2)) * pkin(2) + t26 * t10 - t32 * t9;
t15 = -qJ(3) - pkin(7);
t6 = t10 * t27 + t30;
t5 = t10 * t28 - t29;
t4 = t10 * t29 - t28;
t3 = t10 * t30 + t27;
t2 = t5 * t16 + t6 * t19;
t1 = -t6 * t16 + t5 * t19;
t7 = [-t20 * t15 - t33 * t18 - t22 * t4 - t23 * t3, t21 * t20, t18, t20 * t9, t5, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t5 * qJ(5) - t18 * t15 + t33 * t20 + t31 * t6, t21 * t18, -t20, t18 * t9, t3 (-t4 * t16 + t3 * t19) * r_i_i_C(1) + (-t3 * t16 - t4 * t19) * r_i_i_C(2); 0, t32 * t10 + t34, 0, -t10, t9 * t13 ((t13 * t19 - t14 * t16) * r_i_i_C(1) + (-t13 * t16 - t14 * t19) * r_i_i_C(2)) * t9;];
Ja_transl  = t7;
