% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRP7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:33:37
% EndTime: 2019-02-26 20:33:37
% DurationCPUTime: 0.20s
% Computational Cost: add. (234->53), mult. (388->78), div. (0->0), fcn. (304->7), ass. (0->42)
t216 = sin(qJ(5));
t212 = pkin(9) + qJ(4);
t210 = cos(t212);
t217 = sin(qJ(1));
t209 = sin(t212);
t229 = qJD(1) * t209 + qJD(5);
t219 = cos(qJ(1));
t239 = qJD(4) * t219;
t250 = -t210 * t239 + t229 * t217;
t259 = t250 * t216;
t218 = cos(qJ(5));
t245 = t218 * pkin(5);
t208 = pkin(4) + t245;
t244 = r_i_i_C(3) + qJ(6) + pkin(8);
t251 = t244 * t210;
t258 = (t251 - pkin(3) * sin(pkin(9)) - t208 * t209 - qJ(2)) * qJD(1) - qJD(5) * t245 - qJD(3);
t247 = r_i_i_C(2) * t216;
t226 = r_i_i_C(1) * t218 + t208 - t247;
t257 = -(-t226 * t209 + t251) * qJD(4) - qJD(6) * t209;
t232 = t244 * t209;
t256 = t226 * t210 + t232;
t249 = pkin(5) + r_i_i_C(1);
t225 = r_i_i_C(2) * t218 + t249 * t216;
t248 = pkin(5) * t216;
t243 = t218 * t219;
t242 = qJD(1) * t217;
t211 = qJD(1) * t219;
t241 = qJD(4) * t210;
t240 = qJD(4) * t217;
t238 = qJD(5) * t209;
t237 = qJD(5) * t210;
t235 = t210 * qJD(6);
t230 = qJD(1) + t238;
t227 = t230 * t219;
t223 = qJD(5) * t225;
t221 = qJD(1) * t256;
t206 = -t230 * t218 * t217 + (-t210 * t240 - t229 * t219) * t216;
t220 = -t238 * t248 - t235 + qJD(2) + (t208 * t210 + t232) * qJD(4) + (-pkin(1) - pkin(7) - qJ(3) - t248) * qJD(1);
t207 = t229 * t243 + (-t230 * t216 + t218 * t241) * t217;
t205 = -t216 * t227 - t218 * t250;
t204 = -t218 * t227 + t259;
t1 = [t205 * r_i_i_C(1) + t204 * r_i_i_C(2) + t258 * t217 + t220 * t219, t211, -t242, t219 * t221 + (-t225 * t237 - t257) * t217, -t207 * r_i_i_C(2) + t249 * t206, t209 * t240 - t210 * t211; t207 * r_i_i_C(1) + t206 * r_i_i_C(2) + t220 * t217 - t258 * t219, t242, t211, t217 * t221 + (t210 * t223 + t257) * t219, -t204 * r_i_i_C(1) + t205 * r_i_i_C(2) + (t230 * t243 - t259) * pkin(5), -t209 * t239 - t210 * t242; 0, 0, 0, -t256 * qJD(4) + t209 * t223 + t235 (-t249 * t218 + t247) * t237 + t225 * t209 * qJD(4), t241;];
JaD_transl  = t1;
