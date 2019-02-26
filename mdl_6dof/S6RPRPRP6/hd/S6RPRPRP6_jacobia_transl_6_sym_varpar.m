% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRP6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:46:27
% EndTime: 2019-02-26 20:46:27
% DurationCPUTime: 0.10s
% Computational Cost: add. (115->29), mult. (128->39), div. (0->0), fcn. (144->7), ass. (0->23)
t12 = sin(qJ(5));
t20 = pkin(3) + r_i_i_C(3) + qJ(6) + pkin(8);
t9 = pkin(9) + qJ(3);
t8 = cos(t9);
t18 = t20 * t8;
t7 = sin(t9);
t28 = -(pkin(5) * t12 + qJ(4)) * t7 - cos(pkin(9)) * pkin(2) - pkin(1) - t18;
t26 = pkin(5) + r_i_i_C(1);
t14 = cos(qJ(5));
t25 = pkin(5) * t14 + pkin(4) + pkin(7) + qJ(2);
t15 = cos(qJ(1));
t24 = t12 * t15;
t13 = sin(qJ(1));
t23 = t13 * t12;
t22 = t13 * t14;
t21 = t14 * t15;
t1 = t7 * t21 - t23;
t3 = t7 * t22 + t24;
t17 = r_i_i_C(2) * t14 + t26 * t12 + qJ(4);
t16 = t17 * t8 - t20 * t7;
t4 = -t7 * t23 + t21;
t2 = t7 * t24 + t22;
t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t28 * t13 + t25 * t15, t13, t16 * t15, t15 * t7, -t2 * r_i_i_C(2) + t26 * t1, t15 * t8; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t25 * t13 - t28 * t15, -t15, t16 * t13, t13 * t7, t4 * r_i_i_C(2) + t26 * t3, t13 * t8; 0, 0, t17 * t7 + t18, -t8 (r_i_i_C(2) * t12 - t26 * t14) * t8, t7;];
Ja_transl  = t5;
