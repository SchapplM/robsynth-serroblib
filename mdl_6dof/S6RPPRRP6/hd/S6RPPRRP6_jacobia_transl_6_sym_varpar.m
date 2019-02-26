% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRP6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_jacobia_transl_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:33:08
% EndTime: 2019-02-26 20:33:08
% DurationCPUTime: 0.10s
% Computational Cost: add. (66->25), mult. (146->33), div. (0->0), fcn. (171->6), ass. (0->23)
t10 = cos(qJ(4));
t19 = pkin(8) + r_i_i_C(2);
t14 = t19 * t10;
t7 = sin(qJ(4));
t24 = -pkin(4) * t7 - pkin(1) - qJ(3) + t14;
t15 = r_i_i_C(3) + qJ(6);
t20 = pkin(5) + r_i_i_C(1);
t6 = sin(qJ(5));
t9 = cos(qJ(5));
t23 = t15 * t6 + t20 * t9 + pkin(4);
t8 = sin(qJ(1));
t22 = t8 * t6;
t21 = t8 * t9;
t11 = cos(qJ(1));
t18 = t11 * t6;
t17 = t11 * t9;
t16 = -pkin(7) + qJ(2);
t12 = t23 * t10 + t19 * t7;
t4 = t7 * t17 - t22;
t3 = t7 * t18 + t21;
t2 = t7 * t21 + t18;
t1 = t7 * t22 - t17;
t5 = [-t15 * t1 + t16 * t11 - t20 * t2 + t24 * t8, t8, t11, t12 * t11, t15 * t4 - t20 * t3, t3; -t24 * t11 + t15 * t3 + t16 * t8 + t20 * t4, -t11, t8, t12 * t8, -t20 * t1 + t15 * t2, t1; 0, 0, 0, -t23 * t7 + t14 (t15 * t9 - t20 * t6) * t10, t10 * t6;];
Ja_transl  = t5;
