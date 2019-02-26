% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRP1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:30:10
% EndTime: 2019-02-26 20:30:10
% DurationCPUTime: 0.10s
% Computational Cost: add. (143->31), mult. (110->39), div. (0->0), fcn. (123->9), ass. (0->25)
t24 = r_i_i_C(3) + qJ(6) + pkin(8);
t11 = pkin(10) + qJ(4);
t7 = sin(t11);
t20 = t24 * t7;
t16 = cos(qJ(5));
t6 = pkin(5) * t16 + pkin(4);
t9 = cos(t11);
t28 = t20 + t9 * t6 + cos(pkin(10)) * pkin(3) + pkin(2);
t27 = pkin(5) + r_i_i_C(1);
t12 = qJ(1) + pkin(9);
t8 = sin(t12);
t26 = t16 * t8;
t15 = sin(qJ(5));
t25 = t8 * t15;
t10 = cos(t12);
t23 = t10 * t15;
t22 = t10 * t16;
t19 = pkin(5) * t15 + pkin(7) + qJ(3);
t18 = r_i_i_C(1) * t16 - r_i_i_C(2) * t15 + t6;
t1 = t9 * t25 + t22;
t3 = -t9 * t23 + t26;
t17 = -t18 * t7 + t24 * t9;
t4 = t9 * t22 + t25;
t2 = -t9 * t26 + t23;
t5 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t19 * t10 - t28 * t8, 0, t8, t17 * t10, -t4 * r_i_i_C(2) + t27 * t3, t10 * t7; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t19 * t8 + t28 * t10, 0, -t10, t17 * t8, t2 * r_i_i_C(2) - t27 * t1, t8 * t7; 0, 1, 0, t18 * t9 + t20 (-r_i_i_C(2) * t16 - t27 * t15) * t7, -t9;];
Ja_transl  = t5;
