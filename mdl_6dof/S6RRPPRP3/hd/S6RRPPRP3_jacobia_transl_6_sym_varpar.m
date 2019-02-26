% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRP3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_jacobia_transl_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:08
% EndTime: 2019-02-26 21:26:08
% DurationCPUTime: 0.10s
% Computational Cost: add. (71->27), mult. (137->38), div. (0->0), fcn. (153->6), ass. (0->23)
t11 = cos(qJ(2));
t17 = pkin(2) + pkin(3) + r_i_i_C(3) + qJ(6) + pkin(8);
t15 = t17 * t11;
t10 = cos(qJ(5));
t18 = pkin(5) * t10 + pkin(4) + qJ(3);
t8 = sin(qJ(2));
t25 = -t18 * t8 - pkin(1) - t15;
t24 = pkin(5) + r_i_i_C(1);
t9 = sin(qJ(1));
t22 = t9 * t8;
t12 = cos(qJ(1));
t21 = t12 * t8;
t20 = t9 * t10;
t19 = t10 * t12;
t7 = sin(qJ(5));
t16 = -pkin(5) * t7 + pkin(7) - qJ(4);
t1 = t7 * t22 - t19;
t3 = -t7 * t21 - t20;
t14 = r_i_i_C(1) * t10 - r_i_i_C(2) * t7 + t18;
t13 = t14 * t11 - t17 * t8;
t4 = t8 * t19 - t9 * t7;
t2 = -t12 * t7 - t8 * t20;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t16 * t12 + t25 * t9, t13 * t12, t21, -t9, -t4 * r_i_i_C(2) + t24 * t3, t12 * t11; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t25 * t12 + t16 * t9, t13 * t9, t22, t12, t2 * r_i_i_C(2) - t24 * t1, t9 * t11; 0, t14 * t8 + t15, -t11, 0 (r_i_i_C(2) * t10 + t24 * t7) * t11, t8;];
Ja_transl  = t5;
