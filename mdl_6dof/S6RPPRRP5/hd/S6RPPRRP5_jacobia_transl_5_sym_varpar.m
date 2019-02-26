% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRP5
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
% Datum: 2019-02-26 20:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRP5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_jacobia_transl_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:32:32
% EndTime: 2019-02-26 20:32:32
% DurationCPUTime: 0.09s
% Computational Cost: add. (40->23), mult. (87->32), div. (0->0), fcn. (99->6), ass. (0->21)
t18 = pkin(8) + r_i_i_C(3);
t9 = cos(qJ(4));
t14 = t18 * t9;
t6 = sin(qJ(4));
t21 = -t6 * pkin(4) - pkin(1) - qJ(3) + t14;
t5 = sin(qJ(5));
t7 = sin(qJ(1));
t20 = t7 * t5;
t8 = cos(qJ(5));
t19 = t7 * t8;
t10 = cos(qJ(1));
t17 = t10 * t5;
t16 = t10 * t8;
t15 = -pkin(7) + qJ(2);
t12 = r_i_i_C(1) * t8 - r_i_i_C(2) * t5 + pkin(4);
t11 = t12 * t9 + t18 * t6;
t4 = t6 * t16 - t20;
t3 = -t6 * t17 - t19;
t2 = -t6 * t19 - t17;
t1 = t6 * t20 - t16;
t13 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * t10 + t21 * t7, t7, t10, t11 * t10, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t21 * t10 + t15 * t7, -t10, t7, t11 * t7, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, 0, 0, -t12 * t6 + t14 (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t9, 0;];
Ja_transl  = t13;
