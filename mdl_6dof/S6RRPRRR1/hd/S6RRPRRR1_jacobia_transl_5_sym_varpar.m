% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:54:00
% EndTime: 2019-02-26 21:54:01
% DurationCPUTime: 0.08s
% Computational Cost: add. (129->16), mult. (66->20), div. (0->0), fcn. (68->10), ass. (0->16)
t13 = qJ(2) + pkin(11);
t10 = qJ(4) + t13;
t9 = qJ(5) + t10;
t5 = sin(t9);
t6 = cos(t9);
t21 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
t20 = t21 + pkin(4) * cos(t10);
t26 = t20 + cos(qJ(2)) * pkin(2) + pkin(3) * cos(t13);
t19 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
t16 = t19 - pkin(4) * sin(t10);
t22 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(3) + pkin(7);
t18 = pkin(1) + t26;
t17 = -pkin(3) * sin(t13) - sin(qJ(2)) * pkin(2) + t16;
t15 = cos(qJ(1));
t14 = sin(qJ(1));
t1 = [-t18 * t14 + t22 * t15, t17 * t15, t14, t16 * t15, t19 * t15, 0; t22 * t14 + t18 * t15, t17 * t14, -t15, t16 * t14, t19 * t14, 0; 0, t26, 0, t20, t21, 0;];
Ja_transl  = t1;
