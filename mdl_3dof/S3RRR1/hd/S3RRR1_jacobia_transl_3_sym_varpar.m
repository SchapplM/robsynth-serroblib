% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S3RRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
%
% Output:
% Ja_transl [3x3]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S3RRR1_jacobia_transl_3_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_jacobia_transl_3_sym_varpar: qJ has to be [3x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S3RRR1_jacobia_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_jacobia_transl_3_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:16:20
% EndTime: 2019-02-26 19:16:20
% DurationCPUTime: 0.06s
% Computational Cost: add. (40->8), mult. (18->8), div. (0->0), fcn. (18->6), ass. (0->9)
t7 = qJ(1) + qJ(2);
t6 = qJ(3) + t7;
t3 = sin(t6);
t4 = cos(t6);
t11 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
t10 = t11 + pkin(2) * cos(t7);
t9 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
t8 = -pkin(2) * sin(t7) + t9;
t1 = [-sin(qJ(1)) * pkin(1) + t8, t8, t9; cos(qJ(1)) * pkin(1) + t10, t10, t11; 0, 0, 0;];
Ja_transl  = t1;
