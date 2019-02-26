% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPPR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:59:24
% EndTime: 2019-02-26 19:59:24
% DurationCPUTime: 0.13s
% Computational Cost: add. (160->41), mult. (514->73), div. (0->0), fcn. (488->8), ass. (0->34)
t197 = sin(pkin(10));
t199 = cos(pkin(10));
t204 = cos(qJ(2));
t200 = cos(pkin(6));
t202 = sin(qJ(2));
t218 = t200 * t202;
t192 = t197 * t204 + t199 * t218;
t203 = cos(qJ(3));
t198 = sin(pkin(6));
t201 = sin(qJ(3));
t220 = t198 * t201;
t223 = -t192 * t203 + t199 * t220;
t215 = pkin(3) + pkin(4) - r_i_i_C(2);
t222 = r_i_i_C(1) + qJ(4);
t206 = t222 * t201 + t215 * t203 + pkin(2);
t219 = t198 * t203;
t217 = t200 * t204;
t216 = qJD(2) * t202;
t214 = pkin(8) - r_i_i_C(3) - qJ(5);
t212 = t199 * t217;
t211 = qJD(2) * t198 * t204;
t208 = t197 * t218 - t199 * t204;
t210 = t197 * t220 - t203 * t208;
t209 = t197 * t217 + t199 * t202;
t207 = t200 * t201 + t202 * t219;
t205 = qJD(4) * t201 + (-t215 * t201 + t222 * t203) * qJD(3);
t190 = t208 * qJD(2);
t189 = t209 * qJD(2);
t188 = t192 * qJD(2);
t187 = -qJD(2) * t212 + t197 * t216;
t185 = t207 * qJD(3) + t201 * t211;
t183 = t210 * qJD(3) - t189 * t201;
t181 = -t223 * qJD(3) - t187 * t201;
t1 = [0, qJD(5) * t208 - t214 * t189 + t206 * t190 - t205 * t209, t210 * qJD(4) + t222 * (-t189 * t203 + (t197 * t219 + t201 * t208) * qJD(3)) - t215 * t183, t183, t190, 0; 0, -t192 * qJD(5) - t214 * t187 + t205 * (-t197 * t202 + t212) - t206 * t188, -t223 * qJD(4) + t222 * (-t187 * t203 + (-t192 * t201 - t199 * t219) * qJD(3)) - t215 * t181, t181, -t188, 0; 0 (-qJD(5) * t202 + t205 * t204 + (-t206 * t202 + t214 * t204) * qJD(2)) * t198, t207 * qJD(4) + t222 * (t203 * t211 + (t200 * t203 - t202 * t220) * qJD(3)) - t215 * t185, t185, -t198 * t216, 0;];
JaD_transl  = t1;
