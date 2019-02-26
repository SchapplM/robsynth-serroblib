% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:28:01
% EndTime: 2019-02-26 20:28:01
% DurationCPUTime: 0.13s
% Computational Cost: add. (64->28), mult. (194->42), div. (0->0), fcn. (145->6), ass. (0->18)
t167 = sin(qJ(4));
t169 = cos(qJ(4));
t165 = sin(pkin(9));
t166 = cos(pkin(9));
t175 = r_i_i_C(1) * t166 - r_i_i_C(2) * t165 + pkin(4);
t180 = r_i_i_C(3) + qJ(5);
t171 = -(t180 * t167 + t175 * t169) * qJD(4) + t169 * qJD(5);
t182 = -qJD(3) + t171;
t168 = sin(qJ(1));
t179 = qJD(1) * t168;
t170 = cos(qJ(1));
t164 = qJD(1) * t170;
t178 = qJD(4) * t167;
t176 = qJD(4) * t180;
t174 = t165 * r_i_i_C(1) + t166 * r_i_i_C(2) + pkin(7) - qJ(2);
t173 = -t175 * qJD(4) + qJD(5);
t172 = -t175 * t167 + t180 * t169 - pkin(1) - qJ(3);
t1 = [t170 * qJD(2) + t182 * t168 + (t174 * t168 + t172 * t170) * qJD(1), t164, -t179 (t170 * t176 - t175 * t179) * t169 + (t173 * t170 - t180 * t179) * t167, t169 * t179 + t170 * t178, 0; t168 * qJD(2) - t182 * t170 + (t172 * t168 - t174 * t170) * qJD(1), t179, t164 (t175 * t164 + t168 * t176) * t169 + (t180 * t164 + t173 * t168) * t167, -t169 * t164 + t168 * t178, 0; 0, 0, 0, t171, qJD(4) * t169, 0;];
JaD_transl  = t1;
