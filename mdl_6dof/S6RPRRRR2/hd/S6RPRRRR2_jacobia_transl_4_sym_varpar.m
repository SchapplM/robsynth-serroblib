% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR2_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR2_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:15:35
% EndTime: 2019-02-26 21:15:35
% DurationCPUTime: 0.07s
% Computational Cost: add. (59->12), mult. (43->16), div. (0->0), fcn. (43->8), ass. (0->13)
t9 = qJ(3) + qJ(4);
t5 = sin(t9);
t6 = cos(t9);
t15 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
t18 = t15 + cos(qJ(3)) * pkin(3);
t16 = r_i_i_C(3) + pkin(8) + pkin(7);
t14 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
t13 = pkin(2) + t18;
t12 = -sin(qJ(3)) * pkin(3) + t14;
t8 = qJ(1) + pkin(11);
t4 = cos(t8);
t3 = sin(t8);
t1 = [-sin(qJ(1)) * pkin(1) + t16 * t4 - t13 * t3, 0, t12 * t4, t14 * t4, 0, 0; cos(qJ(1)) * pkin(1) + t16 * t3 + t13 * t4, 0, t12 * t3, t14 * t3, 0, 0; 0, 1, t18, t15, 0, 0;];
Ja_transl  = t1;
