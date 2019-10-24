% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPP3
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:29
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRPP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:31
	% EndTime: 2019-10-24 10:29:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:31
	% EndTime: 2019-10-24 10:29:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:31
	% EndTime: 2019-10-24 10:29:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(2));
	t13 = cos(qJ(2));
	t14 = qJD(2) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, cos(pkin(7)) * t14, 0, 0, 0; 0, sin(pkin(7)) * t14, 0, 0, 0; 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:32
	% EndTime: 2019-10-24 10:29:32
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (31->19), mult. (116->43), div. (0->0), fcn. (87->6), ass. (0->17)
	t151 = pkin(6) + r_i_i_C(3);
	t138 = sin(qJ(3));
	t141 = cos(qJ(2));
	t150 = t138 * t141;
	t140 = cos(qJ(3));
	t149 = t140 * t141;
	t148 = qJD(2) * t141;
	t139 = sin(qJ(2));
	t147 = qJD(3) * t139;
	t146 = r_i_i_C(1) * t138 + r_i_i_C(2) * t140;
	t145 = -r_i_i_C(1) * t140 + r_i_i_C(2) * t138 - pkin(2);
	t144 = t146 * t139;
	t143 = qJD(2) * t144;
	t142 = qJD(3) * t144 + (-t151 * t139 + t145 * t141) * qJD(2);
	t137 = cos(pkin(7));
	t136 = sin(pkin(7));
	t1 = [0, t142 * t137, t137 * t143 + ((-t136 * t138 - t137 * t149) * r_i_i_C(1) + (-t136 * t140 + t137 * t150) * r_i_i_C(2)) * qJD(3), 0, 0; 0, t142 * t136, t136 * t143 + ((-t136 * t149 + t137 * t138) * r_i_i_C(1) + (t136 * t150 + t137 * t140) * r_i_i_C(2)) * qJD(3), 0, 0; 0, -t146 * t141 * qJD(3) + (t145 * t139 + t151 * t141) * qJD(2), (t138 * t147 - t140 * t148) * r_i_i_C(2) + (-t138 * t148 - t140 * t147) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:32
	% EndTime: 2019-10-24 10:29:33
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (89->30), mult. (314->55), div. (0->0), fcn. (266->8), ass. (0->27)
	t213 = sin(qJ(3));
	t215 = cos(qJ(3));
	t209 = sin(pkin(8));
	t211 = cos(pkin(8));
	t224 = t211 * r_i_i_C(1) - t209 * r_i_i_C(2) + pkin(3);
	t232 = r_i_i_C(3) + qJ(4);
	t233 = t224 * t213 - t232 * t215;
	t234 = t233 * qJD(3) - t213 * qJD(4);
	t216 = cos(qJ(2));
	t231 = t213 * t216;
	t230 = t215 * t216;
	t214 = sin(qJ(2));
	t229 = qJD(2) * t214;
	t228 = qJD(2) * t216;
	t226 = t213 * t229;
	t225 = t215 * t229;
	t223 = t209 * r_i_i_C(1) + t211 * r_i_i_C(2) + pkin(6);
	t210 = sin(pkin(7));
	t212 = cos(pkin(7));
	t222 = -t210 * t213 - t212 * t230;
	t221 = -t210 * t230 + t212 * t213;
	t220 = -t232 * t213 - t224 * t215;
	t218 = -pkin(2) + t220;
	t217 = t234 * t214 + (-t223 * t214 + t218 * t216) * qJD(2);
	t207 = t222 * qJD(3) + t212 * t226;
	t205 = t221 * qJD(3) + t210 * t226;
	t1 = [0, t217 * t212, -t222 * qJD(4) - t232 * (t212 * t225 + (-t210 * t215 + t212 * t231) * qJD(3)) + t224 * t207, -t207, 0; 0, t217 * t210, -t221 * qJD(4) - t232 * (t210 * t225 + (t210 * t231 + t212 * t215) * qJD(3)) + t224 * t205, -t205, 0; 0, -t234 * t216 + (t218 * t214 + t223 * t216) * qJD(2), -t233 * t228 + (t220 * qJD(3) + qJD(4) * t215) * t214, t214 * qJD(3) * t215 + t213 * t228, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:32
	% EndTime: 2019-10-24 10:29:33
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (142->53), mult. (497->92), div. (0->0), fcn. (440->8), ass. (0->40)
	t237 = sin(qJ(3));
	t239 = cos(qJ(3));
	t269 = r_i_i_C(2) + qJ(4);
	t252 = t269 * t239;
	t233 = sin(pkin(8));
	t260 = qJD(5) * t233;
	t274 = (-pkin(3) * t237 + t252) * qJD(3) + t237 * qJD(4) + t239 * t260;
	t235 = cos(pkin(8));
	t238 = sin(qJ(2));
	t240 = cos(qJ(2));
	t262 = qJD(3) * t237;
	t255 = t238 * t262;
	t265 = t239 * t240;
	t244 = -t233 * t255 + (t233 * t265 - t238 * t235) * qJD(2);
	t253 = t269 * t237;
	t247 = -pkin(3) * t239 - pkin(2) - t253;
	t259 = t235 * qJD(5);
	t268 = r_i_i_C(3) + qJ(5);
	t271 = pkin(4) + r_i_i_C(1);
	t272 = -t240 * t259 - t274 * t238 + (-t238 * pkin(6) + t247 * t240) * qJD(2) - t268 * t244 + t271 * (t235 * t255 + (-t233 * t238 - t235 * t265) * qJD(2));
	t245 = t268 * t233 + t271 * t235 + pkin(3);
	t267 = t237 * t240;
	t266 = t238 * t239;
	t264 = qJD(2) * t238;
	t263 = qJD(2) * t240;
	t261 = qJD(3) * t239;
	t258 = t237 * t264;
	t257 = t239 * t264;
	t256 = t235 * t263;
	t254 = t240 * t262;
	t234 = sin(pkin(7));
	t236 = cos(pkin(7));
	t251 = -t234 * t237 - t236 * t265;
	t250 = -t234 * t265 + t236 * t237;
	t249 = t234 * t267 + t236 * t239;
	t229 = -t234 * t261 + (t254 + t257) * t236;
	t228 = t251 * qJD(3) + t236 * t258;
	t227 = t249 * qJD(3) + t234 * t257;
	t226 = t250 * qJD(3) + t234 * t258;
	t1 = [0, t272 * t236, (t234 * t239 - t236 * t267) * t260 - t251 * qJD(4) - t269 * t229 + t245 * t228, -t228, -t229 * t233 - t236 * t256; 0, t272 * t234, -t250 * qJD(4) + t245 * t226 - t269 * t227 - t249 * t260, -t226, -t227 * t233 - t234 * t256; 0, -t238 * t259 + t271 * (-t235 * t254 + (t233 * t240 - t235 * t266) * qJD(2)) - t268 * (t233 * t254 + (t233 * t266 + t235 * t240) * qJD(2)) + t274 * t240 + (t240 * pkin(6) + t247 * t238) * qJD(2), (-t237 * t245 + t252) * t263 + (-t237 * t260 + qJD(4) * t239 + (-t239 * t245 - t253) * qJD(3)) * t238, t237 * t263 + t238 * t261, t244;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end