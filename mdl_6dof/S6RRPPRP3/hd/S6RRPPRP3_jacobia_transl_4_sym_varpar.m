% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRP3_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP3_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_jacobia_transl_4_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:08
% EndTime: 2019-02-26 21:26:08
% DurationCPUTime: 0.07s
% Computational Cost: add. (30->10), mult. (55->12), div. (0->0), fcn. (60->4), ass. (0->11)
t1 = sin(qJ(2));
t3 = cos(qJ(2));
t8 = pkin(2) + pkin(3) - r_i_i_C(2);
t9 = r_i_i_C(1) + qJ(3);
t6 = t9 * t1 + t8 * t3;
t10 = pkin(1) + t6;
t7 = pkin(7) - r_i_i_C(3) - qJ(4);
t5 = -t8 * t1 + t9 * t3;
t4 = cos(qJ(1));
t2 = sin(qJ(1));
t11 = [-t10 * t2 + t7 * t4, t5 * t4, t4 * t1, -t2, 0, 0; t10 * t4 + t7 * t2, t5 * t2, t2 * t1, t4, 0, 0; 0, t6, -t3, 0, 0, 0;];
Ja_transl  = t11;