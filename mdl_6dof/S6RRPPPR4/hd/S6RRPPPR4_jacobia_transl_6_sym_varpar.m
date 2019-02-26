% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPPR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:23:45
% EndTime: 2019-02-26 21:23:45
% DurationCPUTime: 0.13s
% Computational Cost: add. (95->36), mult. (230->55), div. (0->0), fcn. (279->8), ass. (0->26)
t11 = sin(qJ(2));
t14 = cos(qJ(2));
t22 = r_i_i_C(3) + pkin(8) - qJ(4) - pkin(2);
t20 = t22 * t14;
t27 = -t11 * qJ(3) - pkin(1) + t20;
t26 = pkin(3) + pkin(7);
t25 = pkin(4) + pkin(5);
t12 = sin(qJ(1));
t24 = t12 * t11;
t15 = cos(qJ(1));
t23 = t15 * t11;
t10 = sin(qJ(6));
t13 = cos(qJ(6));
t19 = t10 * r_i_i_C(1) + t13 * r_i_i_C(2) + qJ(5);
t18 = t13 * r_i_i_C(1) - t10 * r_i_i_C(2) + t25;
t8 = sin(pkin(9));
t9 = cos(pkin(9));
t17 = t18 * t8 - t19 * t9 + qJ(3);
t16 = t22 * t11 + t17 * t14;
t6 = t15 * t9 - t8 * t24;
t5 = t15 * t8 + t9 * t24;
t4 = t12 * t9 + t8 * t23;
t3 = t12 * t8 - t9 * t23;
t2 = t3 * t10 + t4 * t13;
t1 = -t4 * t10 + t3 * t13;
t7 = [t27 * t12 + t26 * t15 + t18 * t6 + t19 * t5, t16 * t15, t23, t15 * t14, t3, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t3 * qJ(5) + t26 * t12 - t27 * t15 + t25 * t4, t16 * t12, t24, t12 * t14, -t5 (t6 * t10 - t5 * t13) * r_i_i_C(1) + (t5 * t10 + t6 * t13) * r_i_i_C(2); 0, t17 * t11 - t20, -t14, t11, t14 * t9 ((t10 * t8 + t13 * t9) * r_i_i_C(1) + (-t10 * t9 + t13 * t8) * r_i_i_C(2)) * t14;];
Ja_transl  = t7;
