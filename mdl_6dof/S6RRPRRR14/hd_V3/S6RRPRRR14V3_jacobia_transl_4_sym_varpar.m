% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR14V3_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14V3_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobia_transl_4_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:04
% EndTime: 2019-04-12 15:12:04
% DurationCPUTime: 0.08s
% Computational Cost: add. (27->16), mult. (72->33), div. (0->0), fcn. (83->6), ass. (0->18)
t6 = sin(qJ(2));
t7 = sin(qJ(1));
t17 = t7 * t6;
t9 = cos(qJ(2));
t16 = t7 * t9;
t10 = cos(qJ(1));
t15 = t10 * t9;
t14 = r_i_i_C(3) + qJ(3);
t13 = t14 * t6;
t5 = sin(qJ(4));
t8 = cos(qJ(4));
t12 = r_i_i_C(1) * t8 - r_i_i_C(2) * t5;
t11 = -t12 * t6 + t14 * t9;
t4 = t8 * t15 + t7 * t5;
t3 = -t5 * t15 + t7 * t8;
t2 = t10 * t5 - t8 * t16;
t1 = t10 * t8 + t5 * t16;
t18 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t14 * t17, t11 * t10, t10 * t6, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t10 * t13, t11 * t7, t17, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0; 0, t12 * t9 + t13, -t9 (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t6, 0, 0;];
Ja_transl  = t18;
