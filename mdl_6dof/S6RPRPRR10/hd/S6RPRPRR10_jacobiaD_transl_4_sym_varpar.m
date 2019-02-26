% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR10_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR10_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:00
% EndTime: 2019-02-26 20:54:00
% DurationCPUTime: 0.13s
% Computational Cost: add. (59->23), mult. (186->39), div. (0->0), fcn. (139->6), ass. (0->19)
t164 = sin(qJ(3));
t166 = cos(qJ(3));
t162 = sin(pkin(10));
t163 = cos(pkin(10));
t170 = r_i_i_C(1) * t163 - r_i_i_C(2) * t162 + pkin(3);
t177 = r_i_i_C(3) + qJ(4);
t181 = -(t177 * t164 + t170 * t166) * qJD(3) + t166 * qJD(4);
t180 = t170 * qJD(3) - qJD(4);
t178 = t170 * t164 - t177 * t166 + qJ(2);
t165 = sin(qJ(1));
t176 = qJD(1) * t165;
t167 = cos(qJ(1));
t175 = qJD(1) * t167;
t174 = qJD(3) * t165;
t173 = qJD(3) * t167;
t171 = qJD(1) * t177;
t169 = qJD(1) * t170;
t168 = qJD(2) + (-t162 * r_i_i_C(1) - t163 * r_i_i_C(2) - pkin(1) - pkin(7)) * qJD(1) - t181;
t1 = [t168 * t167 - t178 * t176, t175 (t167 * t169 + t177 * t174) * t166 + (-t180 * t165 + t167 * t171) * t164, t164 * t174 - t166 * t175, 0, 0; t168 * t165 + t178 * t175, t176 (t165 * t169 - t177 * t173) * t166 + (t165 * t171 + t180 * t167) * t164, -t164 * t173 - t166 * t176, 0, 0; 0, 0, t181, qJD(3) * t166, 0, 0;];
JaD_transl  = t1;
