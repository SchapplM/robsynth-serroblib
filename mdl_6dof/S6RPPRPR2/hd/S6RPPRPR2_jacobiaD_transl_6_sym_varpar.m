% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:26:17
% EndTime: 2019-02-26 20:26:17
% DurationCPUTime: 0.23s
% Computational Cost: add. (301->49), mult. (328->81), div. (0->0), fcn. (256->9), ass. (0->34)
t217 = pkin(10) + qJ(4);
t213 = sin(t217);
t215 = cos(t217);
t234 = pkin(4) + pkin(8) + r_i_i_C(3);
t231 = t234 * t213;
t247 = (-qJ(5) * t215 + t231) * qJD(4) - t213 * qJD(5);
t220 = sin(qJ(6));
t229 = qJD(6) * t213 + qJD(1);
t221 = cos(qJ(6));
t237 = qJD(4) * t221;
t245 = -t215 * t237 + t229 * t220;
t238 = qJD(4) * t220;
t244 = t215 * t238 + t229 * t221;
t243 = pkin(5) + pkin(7) + qJ(3);
t218 = qJ(1) + pkin(9);
t214 = sin(t218);
t241 = qJD(1) * t214;
t216 = cos(t218);
t240 = qJD(1) * t216;
t239 = qJD(4) * t216;
t236 = qJD(6) * t215;
t230 = t234 * t215;
t228 = -qJD(1) * t213 - qJD(6);
t227 = t228 * t220;
t226 = t228 * t221;
t225 = r_i_i_C(1) * t220 + r_i_i_C(2) * t221 + qJ(5);
t224 = -qJ(5) * t213 - cos(pkin(10)) * pkin(3) - pkin(2) - t230;
t223 = qJD(5) + (r_i_i_C(1) * t221 - r_i_i_C(2) * t220) * qJD(6);
t222 = t225 * t215 - t231;
t211 = t214 * t227 + t244 * t216;
t210 = t214 * t226 - t245 * t216;
t209 = -t244 * t214 + t216 * t227;
t208 = t245 * t214 + t216 * t226;
t1 = [t209 * r_i_i_C(1) + t208 * r_i_i_C(2) + t216 * qJD(3) + t247 * t214 + (-cos(qJ(1)) * pkin(1) - t243 * t214 + t224 * t216) * qJD(1), 0, t240 (-t225 * t239 + t234 * t241) * t213 + (-t225 * t241 + (-t234 * qJD(4) + t223) * t216) * t215, -t213 * t241 + t215 * t239, t210 * r_i_i_C(1) - t211 * r_i_i_C(2); t211 * r_i_i_C(1) + t210 * r_i_i_C(2) + t214 * qJD(3) - t247 * t216 + (-sin(qJ(1)) * pkin(1) + t243 * t216 + t224 * t214) * qJD(1), 0, t241, t222 * t240 + (t223 * t215 + (-t225 * t213 - t230) * qJD(4)) * t214, t214 * qJD(4) * t215 + t213 * t240, -t208 * r_i_i_C(1) + t209 * r_i_i_C(2); 0, 0, 0, t222 * qJD(4) + t223 * t213, qJD(4) * t213 (-t213 * t238 + t221 * t236) * r_i_i_C(2) + (t213 * t237 + t220 * t236) * r_i_i_C(1);];
JaD_transl  = t1;
