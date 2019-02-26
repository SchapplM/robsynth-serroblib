% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP8_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP8_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:47:27
% EndTime: 2019-02-26 20:47:27
% DurationCPUTime: 0.18s
% Computational Cost: add. (172->44), mult. (308->72), div. (0->0), fcn. (238->8), ass. (0->34)
t211 = qJ(3) + pkin(9);
t208 = sin(t211);
t209 = cos(t211);
t239 = pkin(8) + r_i_i_C(3);
t224 = t239 * t209 - sin(qJ(3)) * pkin(3);
t213 = sin(qJ(5));
t216 = cos(qJ(5));
t226 = r_i_i_C(1) * t216 - r_i_i_C(2) * t213 + pkin(4);
t248 = (-t226 * t208 + t224) * qJD(3);
t247 = -pkin(4) * t208 - qJ(2) + t224;
t222 = t239 * t208 + cos(qJ(3)) * pkin(3);
t245 = t226 * t209 + t222;
t218 = cos(qJ(1));
t228 = qJD(1) * t208 + qJD(5);
t242 = t228 * t218;
t215 = sin(qJ(1));
t233 = qJD(3) * t209;
t240 = t228 * t215 - t218 * t233;
t235 = -pkin(1) - qJ(4) - pkin(7);
t234 = qJD(1) * t215;
t232 = qJD(3) * t216;
t231 = qJD(5) * t209;
t229 = -qJD(5) * t208 - qJD(1);
t227 = r_i_i_C(1) * t213 + r_i_i_C(2) * t216;
t225 = t229 * t218;
t221 = qJD(5) * t227;
t220 = qJD(2) + (pkin(4) * t209 + t222) * qJD(3);
t219 = qJD(1) * t245;
t210 = qJD(1) * t218;
t207 = t216 * t242 + (t209 * t232 + t229 * t213) * t215;
t206 = t229 * t216 * t215 + (-t215 * t233 - t242) * t213;
t205 = t213 * t225 - t240 * t216;
t204 = t240 * t213 + t216 * t225;
t1 = [t205 * r_i_i_C(1) + t204 * r_i_i_C(2) - t215 * qJD(4) + t220 * t218 + (t247 * t215 + t235 * t218) * qJD(1), t210, t218 * t219 + (-t227 * t231 + t248) * t215, -t234, t206 * r_i_i_C(1) - t207 * r_i_i_C(2), 0; t207 * r_i_i_C(1) + t206 * r_i_i_C(2) + t218 * qJD(4) + t220 * t215 + (t235 * t215 - t247 * t218) * qJD(1), t234, t215 * t219 + (t209 * t221 - t248) * t218, t210, -t204 * r_i_i_C(1) + t205 * r_i_i_C(2), 0; 0, 0, -t245 * qJD(3) + t208 * t221, 0 (t208 * t232 + t213 * t231) * r_i_i_C(2) + (qJD(3) * t208 * t213 - t216 * t231) * r_i_i_C(1), 0;];
JaD_transl  = t1;
