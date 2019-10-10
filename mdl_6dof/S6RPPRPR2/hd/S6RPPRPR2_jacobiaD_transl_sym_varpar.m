% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRPR2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:36
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(9);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (26->10), mult. (32->14), div. (0->0), fcn. (20->6), ass. (0->6)
	t19 = r_i_i_C(3) + qJ(3);
	t18 = -r_i_i_C(1) * cos(pkin(10)) + r_i_i_C(2) * sin(pkin(10)) - pkin(2);
	t15 = qJ(1) + pkin(9);
	t14 = cos(t15);
	t13 = sin(t15);
	t1 = [t14 * qJD(3) + (-cos(qJ(1)) * pkin(1) - t19 * t13 + t18 * t14) * qJD(1), 0, qJD(1) * t14, 0, 0, 0; t13 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t19 * t14 + t18 * t13) * qJD(1), 0, qJD(1) * t13, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (71->22), mult. (74->34), div. (0->0), fcn. (48->7), ass. (0->15)
	t40 = r_i_i_C(3) + pkin(7) + qJ(3);
	t31 = qJ(1) + pkin(9);
	t27 = sin(t31);
	t39 = qJD(1) * t27;
	t29 = cos(t31);
	t38 = qJD(1) * t29;
	t37 = qJD(4) * t27;
	t36 = qJD(4) * t29;
	t30 = pkin(10) + qJ(4);
	t26 = sin(t30);
	t28 = cos(t30);
	t35 = r_i_i_C(1) * t26 + r_i_i_C(2) * t28;
	t34 = -r_i_i_C(1) * t28 + r_i_i_C(2) * t26 - cos(pkin(10)) * pkin(3) - pkin(2);
	t33 = t35 * qJD(4);
	t1 = [t29 * qJD(3) + t35 * t37 + (-cos(qJ(1)) * pkin(1) - t40 * t27 + t34 * t29) * qJD(1), 0, t38, (t26 * t36 + t28 * t39) * r_i_i_C(2) + (t26 * t39 - t28 * t36) * r_i_i_C(1), 0, 0; t27 * qJD(3) - t29 * t33 + (-sin(qJ(1)) * pkin(1) + t40 * t29 + t34 * t27) * qJD(1), 0, t39, (t26 * t37 - t28 * t38) * r_i_i_C(2) + (-t26 * t38 - t28 * t37) * r_i_i_C(1), 0, 0; 0, 0, 0, -t33, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:01
	% EndTime: 2019-10-09 23:36:01
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (148->27), mult. (144->39), div. (0->0), fcn. (100->7), ass. (0->17)
	t156 = pkin(10) + qJ(4);
	t152 = sin(t156);
	t154 = cos(t156);
	t167 = r_i_i_C(3) + qJ(5);
	t169 = pkin(4) - r_i_i_C(2);
	t170 = t169 * t152 - t167 * t154;
	t171 = t170 * qJD(4) - t152 * qJD(5);
	t168 = r_i_i_C(1) + pkin(7) + qJ(3);
	t157 = qJ(1) + pkin(9);
	t153 = sin(t157);
	t166 = qJD(1) * t153;
	t155 = cos(t157);
	t165 = qJD(1) * t155;
	t164 = qJD(4) * t155;
	t162 = -t167 * t152 - t169 * t154;
	t160 = -cos(pkin(10)) * pkin(3) - pkin(2) + t162;
	t1 = [t155 * qJD(3) + t171 * t153 + (-cos(qJ(1)) * pkin(1) - t168 * t153 + t160 * t155) * qJD(1), 0, t165, (-t167 * t164 + t169 * t166) * t152 + (-t167 * t166 + (-t169 * qJD(4) + qJD(5)) * t155) * t154, -t152 * t166 + t154 * t164, 0; t153 * qJD(3) - t171 * t155 + (-sin(qJ(1)) * pkin(1) + t168 * t155 + t160 * t153) * qJD(1), 0, t166, -t170 * t165 + (t162 * qJD(4) + qJD(5) * t154) * t153, t153 * qJD(4) * t154 + t152 * t165, 0; 0, 0, 0, -t171, qJD(4) * t152, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:02
	% EndTime: 2019-10-09 23:36:02
	% DurationCPUTime: 0.26s
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
	t1 = [t209 * r_i_i_C(1) + t208 * r_i_i_C(2) + t216 * qJD(3) + t247 * t214 + (-cos(qJ(1)) * pkin(1) - t243 * t214 + t224 * t216) * qJD(1), 0, t240, (-t225 * t239 + t234 * t241) * t213 + (-t225 * t241 + (-t234 * qJD(4) + t223) * t216) * t215, -t213 * t241 + t215 * t239, t210 * r_i_i_C(1) - t211 * r_i_i_C(2); t211 * r_i_i_C(1) + t210 * r_i_i_C(2) + t214 * qJD(3) - t247 * t216 + (-sin(qJ(1)) * pkin(1) + t243 * t216 + t224 * t214) * qJD(1), 0, t241, t222 * t240 + (t223 * t215 + (-t225 * t213 - t230) * qJD(4)) * t214, t214 * qJD(4) * t215 + t213 * t240, -t208 * r_i_i_C(1) + t209 * r_i_i_C(2); 0, 0, 0, t222 * qJD(4) + t223 * t213, qJD(4) * t213, (-t213 * t238 + t221 * t236) * r_i_i_C(2) + (t213 * t237 + t220 * t236) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end