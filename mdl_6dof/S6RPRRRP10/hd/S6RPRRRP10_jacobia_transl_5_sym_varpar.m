% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP10_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP10_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:04
% EndTime: 2019-02-26 21:13:04
% DurationCPUTime: 0.10s
% Computational Cost: add. (96->30), mult. (126->43), div. (0->0), fcn. (140->8), ass. (0->27)
t17 = cos(qJ(3));
t26 = r_i_i_C(3) + pkin(9) + pkin(8);
t31 = t26 * t17;
t14 = sin(qJ(3));
t12 = qJ(4) + qJ(5);
t10 = sin(t12);
t11 = cos(t12);
t16 = cos(qJ(4));
t9 = pkin(4) * t16 + pkin(3);
t21 = r_i_i_C(1) * t11 - r_i_i_C(2) * t10 + t9;
t30 = t26 * t14 + t21 * t17;
t18 = cos(qJ(1));
t15 = sin(qJ(1));
t25 = t14 * t15;
t5 = -t10 * t25 + t11 * t18;
t6 = t10 * t18 + t11 * t25;
t29 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t24 = t14 * t18;
t7 = t10 * t24 + t11 * t15;
t8 = -t10 * t15 + t11 * t24;
t28 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
t13 = sin(qJ(4));
t27 = pkin(4) * t13;
t23 = pkin(1) + pkin(7) + t27;
t22 = -r_i_i_C(1) * t10 - r_i_i_C(2) * t11;
t20 = t14 * t9 + qJ(2) - t31;
t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) - t23 * t15 + t20 * t18, t15, t30 * t15 (-t13 * t25 + t16 * t18) * pkin(4) + t29, t29, 0; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t20 * t15 + t23 * t18, -t18, -t30 * t18 (t13 * t24 + t15 * t16) * pkin(4) + t28, t28, 0; 0, 0, -t21 * t14 + t31 (t22 - t27) * t17, t22 * t17, 0;];
Ja_transl  = t1;
