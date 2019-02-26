% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:56:55
% EndTime: 2019-02-26 20:56:55
% DurationCPUTime: 0.30s
% Computational Cost: add. (407->60), mult. (650->89), div. (0->0), fcn. (542->8), ass. (0->44)
t221 = sin(qJ(3));
t223 = cos(qJ(3));
t243 = t221 * qJD(6);
t220 = sin(qJ(4));
t244 = qJD(5) * t220;
t240 = pkin(8) - r_i_i_C(3) - qJ(6);
t258 = t240 * t223;
t261 = (-pkin(3) * t221 + t258) * qJD(3) + t223 * t244 - t243;
t222 = cos(qJ(4));
t241 = pkin(4) + pkin(5) + r_i_i_C(1);
t255 = r_i_i_C(2) + qJ(5);
t257 = t241 * t220 - t255 * t222;
t260 = t257 * qJD(4) - t244;
t228 = -t255 * t220 - t241 * t222;
t226 = -pkin(3) + t228;
t225 = t226 * t221 + t258;
t219 = qJ(1) + pkin(9);
t217 = sin(t219);
t254 = t217 * t220;
t218 = cos(t219);
t253 = t218 * t220;
t252 = t222 * t223;
t251 = qJD(1) * t217;
t250 = qJD(1) * t218;
t249 = qJD(1) * t221;
t248 = qJD(3) * t221;
t247 = qJD(3) * t223;
t246 = qJD(4) * t220;
t245 = qJD(4) * t222;
t242 = t222 * qJD(5);
t239 = t222 * t251;
t238 = t217 * t248;
t237 = t217 * t246;
t236 = t218 * t245;
t235 = t240 * t221;
t232 = t218 * t252 + t254;
t230 = t217 * t245 + t220 * t250;
t229 = -pkin(3) * t223 - pkin(2) - t235;
t224 = -qJD(6) * t223 + t260 * t221 + (t226 * t223 - t235) * qJD(3);
t206 = t232 * qJD(1) - t222 * t238 - t223 * t237 - t236;
t205 = -t218 * t246 - t220 * t238 + t230 * t223 - t239;
t204 = t223 * t239 + (t222 * t248 + t223 * t246) * t218 - t230;
t203 = t248 * t253 - t223 * t236 - t237 + (t218 * t222 + t223 * t254) * qJD(1);
t1 = [-t218 * t242 - t255 * t205 - t241 * t206 - t261 * t217 + (-cos(qJ(1)) * pkin(1) - t217 * pkin(7) + t229 * t218) * qJD(1), 0, t224 * t218 - t225 * t251, t232 * qJD(5) + t241 * t203 - t255 * t204, -t203, t217 * t249 - t218 * t247; -t217 * t242 - t255 * t203 - t241 * t204 + t261 * t218 + (-sin(qJ(1)) * pkin(1) + t218 * pkin(7) + t229 * t217) * qJD(1), 0, t224 * t217 + t225 * t250 -(-t217 * t252 + t253) * qJD(5) + t255 * t206 - t241 * t205, t205, -t217 * t247 - t218 * t249; 0, 0, t225 * qJD(3) - t260 * t223 - t243, -t257 * t247 + (t228 * qJD(4) + t242) * t221, t220 * t247 + t221 * t245, -t248;];
JaD_transl  = t1;
