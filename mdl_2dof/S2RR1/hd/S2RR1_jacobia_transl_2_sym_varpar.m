% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S2RR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
%
% Output:
% Ja_transl [3x2]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S2RR1_jacobia_transl_2_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_jacobia_transl_2_sym_varpar: qJ has to be [2x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S2RR1_jacobia_transl_2_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_jacobia_transl_2_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:11:20
% EndTime: 2019-02-26 19:11:20
% DurationCPUTime: 0.06s
% Computational Cost: add. (9->6), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->8)
t7 = pkin(1) + r_i_i_C(3);
t1 = sin(qJ(2));
t3 = cos(qJ(2));
t6 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
t5 = r_i_i_C(1) * t1 + r_i_i_C(2) * t3;
t4 = cos(qJ(1));
t2 = sin(qJ(1));
t8 = [t2 * t7 - t4 * t6, t5 * t2; 0, -t6; t2 * t6 + t4 * t7, t5 * t4;];
Ja_transl  = t8;
