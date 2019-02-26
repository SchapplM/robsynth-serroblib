% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPPR5
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPPR5_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR5_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:38
% EndTime: 2019-02-26 20:00:38
% DurationCPUTime: 0.17s
% Computational Cost: add. (124->35), mult. (404->66), div. (0->0), fcn. (384->8), ass. (0->32)
t220 = sin(pkin(10));
t222 = cos(pkin(10));
t225 = sin(qJ(2));
t223 = cos(pkin(6));
t227 = cos(qJ(2));
t237 = t223 * t227;
t246 = -t220 * t225 + t222 * t237;
t224 = sin(qJ(3));
t226 = cos(qJ(3));
t242 = r_i_i_C(3) + qJ(4);
t244 = pkin(3) - r_i_i_C(2);
t245 = t242 * t224 + t244 * t226 + pkin(2);
t243 = pkin(8) + r_i_i_C(1);
t221 = sin(pkin(6));
t240 = t221 * t224;
t239 = t221 * t226;
t238 = t223 * t225;
t235 = qJD(2) * t221 * t227;
t216 = t220 * t227 + t222 * t238;
t234 = -t216 * t226 + t222 * t240;
t231 = t220 * t238 - t222 * t227;
t233 = t220 * t240 - t226 * t231;
t232 = t220 * t237 + t222 * t225;
t230 = t223 * t224 + t225 * t239;
t229 = qJD(2) * t245;
t228 = qJD(4) * t224 + (-t244 * t224 + t242 * t226) * qJD(3);
t213 = t232 * qJD(2);
t211 = t246 * qJD(2);
t209 = t230 * qJD(3) + t224 * t235;
t207 = t233 * qJD(3) - t213 * t224;
t205 = -t234 * qJD(3) + t211 * t224;
t1 = [0, -t243 * t213 - t228 * t232 + t231 * t229, t233 * qJD(4) + t242 * (-t213 * t226 + (t220 * t239 + t224 * t231) * qJD(3)) - t244 * t207, t207, 0, 0; 0, t243 * t211 - t216 * t229 + t228 * t246, -t234 * qJD(4) + t242 * (t211 * t226 + (-t216 * t224 - t222 * t239) * qJD(3)) - t244 * t205, t205, 0, 0; 0 (t228 * t227 + (-t245 * t225 + t243 * t227) * qJD(2)) * t221, t230 * qJD(4) + t242 * (t226 * t235 + (t223 * t226 - t225 * t240) * qJD(3)) - t244 * t209, t209, 0, 0;];
JaD_transl  = t1;
