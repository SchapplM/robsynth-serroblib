% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRRPR2_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:57
% EndTime: 2019-02-26 19:40:57
% DurationCPUTime: 0.08s
% Computational Cost: add. (41->23), mult. (141->56), div. (0->0), fcn. (164->12), ass. (0->29)
t215 = sin(pkin(11));
t221 = cos(pkin(6));
t237 = t215 * t221;
t216 = sin(pkin(7));
t217 = sin(pkin(6));
t236 = t216 * t217;
t235 = t216 * t221;
t220 = cos(pkin(7));
t234 = t217 * t220;
t218 = cos(pkin(12));
t233 = t218 * t220;
t219 = cos(pkin(11));
t232 = t219 * t221;
t224 = cos(qJ(4));
t231 = qJD(3) * t224;
t214 = sin(pkin(12));
t210 = -t215 * t214 + t218 * t232;
t230 = t210 * t220 - t219 * t236;
t212 = -t219 * t214 - t218 * t237;
t229 = t212 * t220 + t215 * t236;
t211 = t214 * t232 + t215 * t218;
t223 = sin(qJ(3));
t225 = cos(qJ(3));
t228 = t211 * t225 + t230 * t223;
t213 = -t214 * t237 + t219 * t218;
t227 = t213 * t225 + t229 * t223;
t226 = t223 * t235 + (t214 * t225 + t223 * t233) * t217;
t222 = sin(qJ(4));
t1 = [0, 0, 0, t227 * qJD(3), 0 (-t227 * t222 + (-t212 * t216 + t215 * t234) * t224) * qJD(4) + (-t213 * t223 + t229 * t225) * t231; 0, 0, 0, t228 * qJD(3), 0 (-t228 * t222 + (-t210 * t216 - t219 * t234) * t224) * qJD(4) + (-t211 * t223 + t230 * t225) * t231; 0, 0, 0, t226 * qJD(3), 0 (-t226 * t222 + (-t218 * t236 + t221 * t220) * t224) * qJD(4) + (t225 * t235 + (-t214 * t223 + t225 * t233) * t217) * t231;];
JgD_rot  = t1;
