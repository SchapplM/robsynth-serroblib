% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRP5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:27:26
% EndTime: 2019-02-26 21:27:26
% DurationCPUTime: 0.10s
% Computational Cost: add. (79->25), mult. (114->38), div. (0->0), fcn. (128->8), ass. (0->21)
t11 = sin(qJ(2));
t13 = cos(qJ(2));
t19 = pkin(2) + r_i_i_C(3) + pkin(8) + qJ(4);
t17 = t19 * t13;
t18 = pkin(4) * sin(pkin(9)) + qJ(3);
t24 = -t18 * t11 - pkin(1) - t17;
t22 = pkin(7) + cos(pkin(9)) * pkin(4) + pkin(3);
t12 = sin(qJ(1));
t21 = t12 * t11;
t14 = cos(qJ(1));
t20 = t14 * t11;
t8 = pkin(9) + qJ(5);
t6 = sin(t8);
t7 = cos(t8);
t16 = r_i_i_C(1) * t6 + r_i_i_C(2) * t7 + t18;
t15 = -t19 * t11 + t16 * t13;
t4 = t14 * t7 - t6 * t21;
t3 = t14 * t6 + t7 * t21;
t2 = t12 * t7 + t6 * t20;
t1 = -t12 * t6 + t7 * t20;
t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t24 * t12 + t22 * t14, t15 * t14, t20, t14 * t13, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t22 * t12 - t24 * t14, t15 * t12, t21, t12 * t13, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, t16 * t11 + t17, -t13, t11 (-r_i_i_C(1) * t7 + r_i_i_C(2) * t6) * t13, 0;];
Ja_transl  = t5;
