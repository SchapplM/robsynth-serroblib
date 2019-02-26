% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR14_jacobia_transl_3_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobia_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_transl_3_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:55:54
% EndTime: 2019-02-26 22:55:54
% DurationCPUTime: 0.13s
% Computational Cost: add. (75->32), mult. (199->56), div. (0->0), fcn. (251->10), ass. (0->26)
t26 = r_i_i_C(3) + qJ(3);
t11 = cos(pkin(14));
t12 = cos(pkin(7));
t8 = sin(pkin(14));
t9 = sin(pkin(7));
t29 = (r_i_i_C(1) * t8 + r_i_i_C(2) * t11) * t12 - t26 * t9;
t10 = sin(pkin(6));
t14 = sin(qJ(1));
t28 = t10 * t14;
t16 = cos(qJ(1));
t27 = t10 * t16;
t25 = cos(pkin(6));
t24 = t14 * t25;
t23 = t16 * t25;
t21 = r_i_i_C(1) * t11 - r_i_i_C(2) * t8 + pkin(2);
t13 = sin(qJ(2));
t15 = cos(qJ(2));
t4 = -t16 * t13 - t15 * t24;
t20 = t12 * t4 + t9 * t28;
t1 = t12 * t28 - t4 * t9;
t2 = t13 * t14 - t15 * t23;
t19 = t12 * t2 + t9 * t27;
t18 = t12 * t27 - t2 * t9;
t5 = -t13 * t24 + t15 * t16;
t3 = -t13 * t23 - t14 * t15;
t6 = [(t3 * t11 + t19 * t8) * r_i_i_C(1) + (t19 * t11 - t3 * t8) * r_i_i_C(2) + t3 * pkin(2) - t14 * pkin(1) + pkin(10) * t27 + t26 * t18, t21 * t4 - t29 * t5, t1, 0, 0, 0; (t11 * t5 + t20 * t8) * r_i_i_C(1) + (t20 * t11 - t5 * t8) * r_i_i_C(2) + t5 * pkin(2) + t16 * pkin(1) + pkin(10) * t28 + t26 * t1, -t21 * t2 + t29 * t3, -t18, 0, 0, 0; 0 (-t13 * t29 + t21 * t15) * t10, -t10 * t15 * t9 + t25 * t12, 0, 0, 0;];
Ja_transl  = t6;
