% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% Ja_transl [3x5]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-12 14:37
% Revision: aab8d7cd0cba739f5e0ec8d53b8419901d1154b0 (2019-06-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRRR1_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobia_transl_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR1_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobia_transl_4_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-12 14:37:33
% EndTime: 2019-06-12 14:37:34
% DurationCPUTime: 0.08s
% Computational Cost: add. (24->17), mult. (63->33), div. (0->0), fcn. (73->6), ass. (0->16)
t6 = sin(qJ(3));
t15 = t6 * r_i_i_C(3);
t7 = sin(qJ(1));
t9 = cos(qJ(3));
t14 = t7 * t9;
t10 = cos(qJ(1));
t13 = t10 * t9;
t5 = sin(qJ(4));
t8 = cos(qJ(4));
t12 = r_i_i_C(1) * t8 - r_i_i_C(2) * t5;
t11 = r_i_i_C(3) * t9 - t12 * t6;
t4 = t8 * t13 + t7 * t5;
t3 = -t5 * t13 + t7 * t8;
t2 = t10 * t5 - t8 * t14;
t1 = t10 * t8 + t5 * t14;
t16 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t10 * qJ(2) - t7 * t15, t7, t11 * t10, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t7 * qJ(2) + t10 * t15, -t10, t11 * t7, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0; 0, 0, t12 * t9 + t15, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t6, 0;];
Ja_transl  = t16;
