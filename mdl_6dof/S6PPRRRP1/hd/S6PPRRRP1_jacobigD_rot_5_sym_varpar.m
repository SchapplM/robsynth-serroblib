% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRRRP1_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_jacobigD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:41:34
% EndTime: 2019-02-26 19:41:34
% DurationCPUTime: 0.08s
% Computational Cost: add. (41->23), mult. (141->56), div. (0->0), fcn. (164->12), ass. (0->29)
t212 = sin(pkin(11));
t218 = cos(pkin(6));
t234 = t212 * t218;
t213 = sin(pkin(7));
t214 = sin(pkin(6));
t233 = t213 * t214;
t232 = t213 * t218;
t217 = cos(pkin(7));
t231 = t214 * t217;
t215 = cos(pkin(12));
t230 = t215 * t217;
t216 = cos(pkin(11));
t229 = t216 * t218;
t219 = sin(qJ(4));
t228 = qJD(3) * t219;
t211 = sin(pkin(12));
t207 = -t212 * t211 + t215 * t229;
t227 = t207 * t217 - t216 * t233;
t209 = -t216 * t211 - t215 * t234;
t226 = t209 * t217 + t212 * t233;
t208 = t211 * t229 + t212 * t215;
t220 = sin(qJ(3));
t222 = cos(qJ(3));
t225 = t208 * t222 + t227 * t220;
t210 = -t211 * t234 + t216 * t215;
t224 = t210 * t222 + t226 * t220;
t223 = t220 * t232 + (t211 * t222 + t220 * t230) * t214;
t221 = cos(qJ(4));
t1 = [0, 0, 0, t224 * qJD(3) (t224 * t221 + (-t209 * t213 + t212 * t231) * t219) * qJD(4) + (-t210 * t220 + t226 * t222) * t228, 0; 0, 0, 0, t225 * qJD(3) (t225 * t221 + (-t207 * t213 - t216 * t231) * t219) * qJD(4) + (-t208 * t220 + t227 * t222) * t228, 0; 0, 0, 0, t223 * qJD(3) (t223 * t221 + (-t215 * t233 + t218 * t217) * t219) * qJD(4) + (t222 * t232 + (-t211 * t220 + t222 * t230) * t214) * t228, 0;];
JgD_rot  = t1;
