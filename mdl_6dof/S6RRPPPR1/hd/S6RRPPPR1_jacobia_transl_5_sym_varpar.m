% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPPR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:58
% EndTime: 2019-02-26 21:21:58
% DurationCPUTime: 0.10s
% Computational Cost: add. (94->25), mult. (111->30), div. (0->0), fcn. (129->8), ass. (0->22)
t20 = r_i_i_C(2) + qJ(4);
t9 = qJ(2) + pkin(9);
t6 = sin(t9);
t28 = t20 * t6 + cos(qJ(2)) * pkin(2);
t7 = cos(t9);
t27 = t7 * pkin(3) + pkin(1) + t28;
t10 = sin(pkin(10));
t11 = cos(pkin(10));
t19 = r_i_i_C(3) + qJ(5);
t25 = pkin(4) + r_i_i_C(1);
t26 = t19 * t10 + t25 * t11 + pkin(3);
t15 = cos(qJ(1));
t24 = t10 * t15;
t23 = t11 * t15;
t14 = sin(qJ(1));
t22 = t14 * t10;
t21 = t14 * t11;
t16 = -sin(qJ(2)) * pkin(2) + t20 * t7 - t26 * t6;
t12 = -qJ(3) - pkin(7);
t3 = t7 * t24 - t21;
t1 = t7 * t22 + t23;
t2 = [-t12 * t15 + t25 * (-t7 * t21 + t24) - t19 * t1 - t27 * t14, t16 * t15, t14, t15 * t6, t3, 0; -t14 * t12 + t25 * (t7 * t23 + t22) + t19 * t3 + t27 * t15, t16 * t14, -t15, t14 * t6, t1, 0; 0, t26 * t7 + t28, 0, -t7, t6 * t10, 0;];
Ja_transl  = t2;
