% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRP9_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP9_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:48:01
% EndTime: 2019-02-26 20:48:01
% DurationCPUTime: 0.09s
% Computational Cost: add. (75->27), mult. (96->38), div. (0->0), fcn. (109->8), ass. (0->21)
t13 = cos(qJ(3));
t20 = r_i_i_C(3) + pkin(8) + qJ(4);
t22 = t20 * t13;
t11 = sin(qJ(3));
t5 = cos(pkin(9)) * pkin(4) + pkin(3);
t8 = pkin(9) + qJ(5);
t6 = sin(t8);
t7 = cos(t8);
t16 = r_i_i_C(1) * t7 - r_i_i_C(2) * t6 + t5;
t21 = t20 * t11 + t16 * t13;
t12 = sin(qJ(1));
t19 = t11 * t12;
t14 = cos(qJ(1));
t18 = t11 * t14;
t17 = pkin(4) * sin(pkin(9)) + pkin(1) + pkin(7);
t15 = t11 * t5 + qJ(2) - t22;
t4 = -t12 * t6 + t7 * t18;
t3 = t12 * t7 + t6 * t18;
t2 = t14 * t6 + t7 * t19;
t1 = t14 * t7 - t6 * t19;
t9 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t17 * t12 + t15 * t14, t12, t21 * t12, -t12 * t13, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * t12 + t17 * t14, -t14, -t21 * t14, t14 * t13, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, 0, -t16 * t11 + t22, t11 (-r_i_i_C(1) * t6 - r_i_i_C(2) * t7) * t13, 0;];
Ja_transl  = t9;
