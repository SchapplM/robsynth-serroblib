% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRP3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:44:41
% EndTime: 2019-02-26 20:44:41
% DurationCPUTime: 0.10s
% Computational Cost: add. (114->27), mult. (96->40), div. (0->0), fcn. (107->10), ass. (0->22)
t15 = cos(qJ(3));
t14 = sin(qJ(3));
t21 = r_i_i_C(3) + pkin(8) + qJ(4);
t18 = t21 * t14;
t5 = cos(pkin(10)) * pkin(4) + pkin(3);
t24 = t15 * t5 + pkin(2) + t18;
t11 = qJ(1) + pkin(9);
t7 = sin(t11);
t23 = t15 * t7;
t9 = cos(t11);
t22 = t15 * t9;
t19 = pkin(4) * sin(pkin(10)) + pkin(7);
t10 = pkin(10) + qJ(5);
t6 = sin(t10);
t8 = cos(t10);
t17 = r_i_i_C(1) * t8 - r_i_i_C(2) * t6 + t5;
t16 = -t17 * t14 + t21 * t15;
t4 = t8 * t22 + t7 * t6;
t3 = -t6 * t22 + t7 * t8;
t2 = -t8 * t23 + t9 * t6;
t1 = t6 * t23 + t9 * t8;
t12 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t19 * t9 - t24 * t7, 0, t16 * t9, t9 * t14, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t19 * t7 + t24 * t9, 0, t16 * t7, t7 * t14, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, 1, t17 * t15 + t18, -t15 (-r_i_i_C(1) * t6 - r_i_i_C(2) * t8) * t14, 0;];
Ja_transl  = t12;
