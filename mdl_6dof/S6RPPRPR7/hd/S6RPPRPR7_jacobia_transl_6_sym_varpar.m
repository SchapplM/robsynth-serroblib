% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRPR7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:29:09
% EndTime: 2019-02-26 20:29:09
% DurationCPUTime: 0.10s
% Computational Cost: add. (113->30), mult. (100->37), div. (0->0), fcn. (115->9), ass. (0->24)
t21 = r_i_i_C(3) + pkin(8) + qJ(5);
t11 = pkin(9) + qJ(4);
t9 = cos(t11);
t27 = t21 * t9;
t5 = cos(pkin(10)) * pkin(5) + pkin(4);
t10 = pkin(10) + qJ(6);
t6 = sin(t10);
t8 = cos(t10);
t19 = r_i_i_C(1) * t8 - r_i_i_C(2) * t6 + t5;
t7 = sin(t11);
t26 = t19 * t9 + t21 * t7;
t16 = sin(qJ(1));
t25 = t16 * t6;
t24 = t16 * t8;
t17 = cos(qJ(1));
t23 = t17 * t6;
t22 = t17 * t8;
t20 = pkin(5) * sin(pkin(10)) + pkin(1) + pkin(7) + qJ(3);
t18 = pkin(3) * sin(pkin(9)) - t27 + t7 * t5 + qJ(2);
t4 = t7 * t22 - t25;
t3 = t7 * t23 + t24;
t2 = t7 * t24 + t23;
t1 = -t7 * t25 + t22;
t12 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t20 * t16 + t18 * t17, t16, t17, t26 * t16, -t16 * t9, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t18 * t16 + t20 * t17, -t17, t16, -t26 * t17, t17 * t9, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 0, 0, -t19 * t7 + t27, t7 (-r_i_i_C(1) * t6 - r_i_i_C(2) * t8) * t9;];
Ja_transl  = t12;
