% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPPRR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:17
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PPPRR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPPRR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPPRR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (13->11), mult. (48->28), div. (0->0), fcn. (46->8), ass. (0->14)
	t71 = sin(pkin(8));
	t76 = sin(qJ(4));
	t80 = t71 * t76;
	t77 = cos(qJ(4));
	t79 = t71 * t77;
	t73 = cos(pkin(9));
	t74 = cos(pkin(8));
	t78 = t73 * t74;
	t75 = cos(pkin(7));
	t72 = sin(pkin(7));
	t70 = sin(pkin(9));
	t69 = t72 * t70 + t75 * t78;
	t68 = -t75 * t70 + t72 * t78;
	t1 = [0, 0, 0, ((-t69 * t77 - t75 * t80) * r_i_i_C(1) + (t69 * t76 - t75 * t79) * r_i_i_C(2)) * qJD(4), 0; 0, 0, 0, ((-t68 * t77 - t72 * t80) * r_i_i_C(1) + (t68 * t76 - t72 * t79) * r_i_i_C(2)) * qJD(4), 0; 0, 0, 0, ((-t73 * t79 + t74 * t76) * r_i_i_C(1) + (t73 * t80 + t74 * t77) * r_i_i_C(2)) * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (92->34), mult. (296->70), div. (0->0), fcn. (307->10), ass. (0->34)
	t248 = -pkin(6) - r_i_i_C(3);
	t226 = sin(pkin(9));
	t227 = sin(pkin(8));
	t247 = t226 * t227;
	t233 = sin(qJ(4));
	t246 = t227 * t233;
	t235 = cos(qJ(4));
	t245 = t227 * t235;
	t228 = sin(pkin(7));
	t230 = cos(pkin(8));
	t244 = t228 * t230;
	t231 = cos(pkin(7));
	t243 = t231 * t226;
	t229 = cos(pkin(9));
	t242 = t231 * t229;
	t241 = qJD(4) * t233;
	t240 = qJD(4) * t245;
	t232 = sin(qJ(5));
	t234 = cos(qJ(5));
	t239 = t232 * r_i_i_C(1) + t234 * r_i_i_C(2);
	t219 = t229 * t244 - t243;
	t213 = t219 * t235 + t228 * t246;
	t221 = t228 * t226 + t230 * t242;
	t215 = t221 * t235 + t231 * t246;
	t223 = t229 * t245 - t230 * t233;
	t238 = t229 * t246 + t230 * t235;
	t237 = qJD(5) * t239;
	t236 = qJD(4) * (t234 * r_i_i_C(1) - t232 * r_i_i_C(2) + pkin(4));
	t220 = -t228 * t229 + t230 * t243;
	t218 = t226 * t244 + t242;
	t216 = t238 * qJD(4);
	t210 = t221 * t241 - t231 * t240;
	t208 = t219 * t241 - t228 * t240;
	t1 = [0, 0, 0, t248 * t210 - (-t221 * t233 + t231 * t245) * t237 - t215 * t236, t239 * t210 + ((-t215 * t234 - t220 * t232) * r_i_i_C(1) + (t215 * t232 - t220 * t234) * r_i_i_C(2)) * qJD(5); 0, 0, 0, t248 * t208 - (-t219 * t233 + t228 * t245) * t237 - t213 * t236, t239 * t208 + ((-t213 * t234 - t218 * t232) * r_i_i_C(1) + (t213 * t232 - t218 * t234) * r_i_i_C(2)) * qJD(5); 0, 0, 0, t248 * t216 - t223 * t236 + t238 * t237, t239 * t216 + ((-t223 * t234 - t232 * t247) * r_i_i_C(1) + (t223 * t232 - t234 * t247) * r_i_i_C(2)) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end