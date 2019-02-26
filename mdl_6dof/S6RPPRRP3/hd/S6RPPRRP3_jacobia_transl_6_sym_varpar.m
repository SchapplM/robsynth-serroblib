% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRP3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:31:17
% EndTime: 2019-02-26 20:31:17
% DurationCPUTime: 0.11s
% Computational Cost: add. (130->29), mult. (146->37), div. (0->0), fcn. (169->8), ass. (0->22)
t10 = sin(qJ(4));
t12 = cos(qJ(4));
t18 = pkin(8) + r_i_i_C(2);
t11 = cos(qJ(5));
t15 = r_i_i_C(3) + qJ(6);
t19 = pkin(5) + r_i_i_C(1);
t9 = sin(qJ(5));
t21 = t19 * t11 + t15 * t9 + pkin(4);
t24 = t18 * t10 + t21 * t12;
t22 = t18 * t12;
t20 = pkin(2) + pkin(7);
t17 = t10 * t9;
t16 = t10 * t11;
t14 = t10 * pkin(4) + qJ(3) - t22;
t8 = qJ(1) + pkin(9);
t7 = cos(t8);
t6 = sin(t8);
t4 = t7 * t16 - t6 * t9;
t3 = t11 * t6 + t7 * t17;
t2 = t6 * t16 + t7 * t9;
t1 = -t7 * t11 + t6 * t17;
t5 = [-sin(qJ(1)) * pkin(1) - t20 * t6 + t19 * t4 + t15 * t3 + t14 * t7, 0, t6, t24 * t6, -t19 * t1 + t15 * t2, t1; cos(qJ(1)) * pkin(1) + t20 * t7 + t19 * t2 + t15 * t1 + t14 * t6, 0, -t7, -t24 * t7, -t15 * t4 + t19 * t3, -t3; 0, 1, 0, -t10 * t21 + t22 (t15 * t11 - t19 * t9) * t12, t12 * t9;];
Ja_transl  = t5;
