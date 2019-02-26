% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_transl = S6RPPRRP5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_jacobia_transl_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:32:32
% EndTime: 2019-02-26 20:32:32
% DurationCPUTime: 0.10s
% Computational Cost: add. (59->28), mult. (110->36), div. (0->0), fcn. (125->6), ass. (0->23)
t11 = cos(qJ(4));
t21 = r_i_i_C(3) + qJ(6) + pkin(8);
t17 = t21 * t11;
t10 = cos(qJ(5));
t5 = pkin(5) * t10 + pkin(4);
t8 = sin(qJ(4));
t24 = -t8 * t5 - pkin(1) - qJ(3) + t17;
t23 = pkin(5) + r_i_i_C(1);
t7 = sin(qJ(5));
t9 = sin(qJ(1));
t22 = t9 * t7;
t12 = cos(qJ(1));
t20 = t12 * t7;
t19 = t9 * t10;
t18 = t10 * t12;
t15 = -pkin(5) * t7 - pkin(7) + qJ(2);
t14 = r_i_i_C(1) * t10 - r_i_i_C(2) * t7 + t5;
t1 = t8 * t22 - t18;
t3 = -t8 * t20 - t19;
t13 = t14 * t11 + t21 * t8;
t4 = t8 * t18 - t22;
t2 = -t8 * t19 - t20;
t6 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * t12 + t24 * t9, t9, t12, t13 * t12, -t4 * r_i_i_C(2) + t23 * t3, -t12 * t11; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t24 * t12 + t15 * t9, -t12, t9, t13 * t9, t2 * r_i_i_C(2) - t23 * t1, -t9 * t11; 0, 0, 0, -t14 * t8 + t17 (-r_i_i_C(2) * t10 - t23 * t7) * t11, t8;];
Ja_transl  = t6;
