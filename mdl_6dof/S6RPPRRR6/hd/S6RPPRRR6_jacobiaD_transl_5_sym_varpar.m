% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:37:36
% EndTime: 2019-02-26 20:37:36
% DurationCPUTime: 0.16s
% Computational Cost: add. (92->44), mult. (286->73), div. (0->0), fcn. (223->6), ass. (0->33)
t206 = cos(qJ(4));
t203 = sin(qJ(4));
t228 = pkin(8) + r_i_i_C(3);
t229 = t228 * t203;
t233 = (pkin(4) * t206 + t229) * qJD(4) + qJD(3);
t202 = sin(qJ(5));
t205 = cos(qJ(5));
t212 = r_i_i_C(1) * t205 - r_i_i_C(2) * t202 + pkin(4);
t231 = t212 * t206 + t229;
t226 = pkin(7) - qJ(2);
t207 = cos(qJ(1));
t225 = t205 * t207;
t204 = sin(qJ(1));
t224 = qJD(1) * t204;
t201 = qJD(1) * t207;
t223 = qJD(4) * t203;
t222 = qJD(4) * t206;
t221 = qJD(4) * t207;
t220 = qJD(5) * t203;
t219 = qJD(5) * t206;
t216 = t228 * t206;
t215 = qJD(1) + t220;
t214 = qJD(1) * t203 + qJD(5);
t213 = r_i_i_C(1) * t202 + r_i_i_C(2) * t205;
t211 = t215 * t202;
t210 = t213 * qJD(5);
t209 = -pkin(4) * t203 - pkin(1) - qJ(3) + t216;
t208 = t214 * t204 - t206 * t221;
t200 = -t214 * t225 + (-t205 * t222 + t211) * t204;
t199 = t215 * t205 * t204 + (t204 * t222 + t214 * t207) * t202;
t198 = t208 * t205 + t207 * t211;
t197 = t208 * t202 - t215 * t225;
t1 = [t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t207 * qJD(2) - t233 * t204 + (t226 * t204 + t209 * t207) * qJD(1), t201, -t224 (-t212 * t221 - t228 * t224) * t203 + (-t212 * t224 + (t228 * qJD(4) - t210) * t207) * t206, t197 * r_i_i_C(1) + t198 * r_i_i_C(2), 0; -t198 * r_i_i_C(1) + t197 * r_i_i_C(2) + t204 * qJD(2) + t233 * t207 + (t209 * t204 - t226 * t207) * qJD(1), t224, t201, t231 * t201 + (-t206 * t210 + (-t212 * t203 + t216) * qJD(4)) * t204, -t199 * r_i_i_C(1) + t200 * r_i_i_C(2), 0; 0, 0, 0, -t231 * qJD(4) + t213 * t220 (t202 * t219 + t205 * t223) * r_i_i_C(2) + (t202 * t223 - t205 * t219) * r_i_i_C(1), 0;];
JaD_transl  = t1;
