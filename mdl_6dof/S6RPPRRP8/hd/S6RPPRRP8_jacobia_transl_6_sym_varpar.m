% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRP8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:34:12
% EndTime: 2019-02-26 20:34:12
% DurationCPUTime: 0.11s
% Computational Cost: add. (118->28), mult. (148->34), div. (0->0), fcn. (173->7), ass. (0->24)
t23 = pkin(8) + r_i_i_C(2);
t11 = sin(qJ(5));
t13 = cos(qJ(5));
t17 = r_i_i_C(3) + qJ(6);
t24 = pkin(5) + r_i_i_C(1);
t25 = t17 * t11 + t24 * t13 + pkin(4);
t8 = pkin(9) + qJ(4);
t6 = sin(t8);
t7 = cos(t8);
t28 = t23 * t6 + t25 * t7;
t26 = t23 * t7;
t22 = pkin(1) + pkin(7) + qJ(3);
t14 = cos(qJ(1));
t21 = t11 * t14;
t12 = sin(qJ(1));
t20 = t12 * t11;
t19 = t12 * t13;
t18 = t14 * t13;
t15 = pkin(3) * sin(pkin(9)) + pkin(4) * t6 - t26 + qJ(2);
t4 = t6 * t18 - t20;
t3 = t6 * t21 + t19;
t2 = t6 * t19 + t21;
t1 = t6 * t20 - t18;
t5 = [-t22 * t12 + t15 * t14 + t17 * t3 + t24 * t4, t12, t14, t28 * t12, -t24 * t1 + t17 * t2, t1; t17 * t1 + t15 * t12 + t22 * t14 + t24 * t2, -t14, t12, -t28 * t14, -t17 * t4 + t24 * t3, -t3; 0, 0, 0, -t25 * t6 + t26 (-t24 * t11 + t17 * t13) * t7, t7 * t11;];
Ja_transl  = t5;
