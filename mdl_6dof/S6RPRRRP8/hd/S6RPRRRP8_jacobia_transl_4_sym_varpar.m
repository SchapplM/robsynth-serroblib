% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRP8
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
% Datum: 2019-02-26 21:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP8_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP8_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_jacobia_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:11:54
% EndTime: 2019-02-26 21:11:54
% DurationCPUTime: 0.07s
% Computational Cost: add. (40->14), mult. (45->16), div. (0->0), fcn. (47->6), ass. (0->15)
t17 = cos(qJ(3)) * pkin(3);
t5 = qJ(3) + qJ(4);
t4 = cos(t5);
t16 = r_i_i_C(1) * t4;
t3 = sin(t5);
t15 = r_i_i_C(2) * t3;
t14 = pkin(1) + r_i_i_C(3) + pkin(8) + pkin(7);
t13 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
t12 = -sin(qJ(3)) * pkin(3) + t13;
t11 = qJ(2) - t12;
t9 = cos(qJ(1));
t7 = sin(qJ(1));
t2 = t9 * t15;
t1 = t7 * t16;
t6 = [t11 * t9 - t14 * t7, t7, t1 + (-t15 + t17) * t7, -t7 * t15 + t1, 0, 0; t11 * t7 + t14 * t9, -t9, t2 + (-t16 - t17) * t9, -t9 * t16 + t2, 0, 0; 0, 0, t12, t13, 0, 0;];
Ja_transl  = t6;
