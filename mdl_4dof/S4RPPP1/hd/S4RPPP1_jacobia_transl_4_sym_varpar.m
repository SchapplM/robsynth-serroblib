% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RPPP1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
%
% Output:
% Ja_transl [3x4]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4RPPP1_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobia_transl_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPPP1_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobia_transl_4_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:29:58
% EndTime: 2019-02-26 19:29:58
% DurationCPUTime: 0.06s
% Computational Cost: add. (32->16), mult. (70->22), div. (0->0), fcn. (93->6), ass. (0->17)
t5 = sin(pkin(6));
t9 = sin(qJ(1));
t16 = t9 * t5;
t7 = cos(pkin(6));
t15 = t9 * t7;
t10 = cos(qJ(1));
t8 = cos(pkin(4));
t14 = t10 * t8;
t13 = r_i_i_C(2) + qJ(3);
t12 = pkin(2) + r_i_i_C(3) + qJ(4);
t6 = sin(pkin(4));
t11 = (pkin(3) + r_i_i_C(1) + qJ(2)) * t6;
t4 = t10 * t7 - t8 * t16;
t3 = t10 * t5 + t8 * t15;
t2 = t5 * t14 + t15;
t1 = -t7 * t14 + t16;
t17 = [-t9 * pkin(1) - t13 * t1 + t10 * t11 - t12 * t2, t9 * t6, t3, t4; t10 * pkin(1) + t9 * t11 + t12 * t4 + t13 * t3, -t10 * t6, t1, t2; 0, t8, -t6 * t7, t6 * t5;];
Ja_transl  = t17;
