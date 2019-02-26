% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRP10_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:01
% EndTime: 2019-02-26 22:45:01
% DurationCPUTime: 0.10s
% Computational Cost: add. (38->16), mult. (130->40), div. (0->0), fcn. (134->8), ass. (0->24)
t220 = sin(pkin(6));
t225 = cos(qJ(3));
t239 = t220 * t225;
t227 = cos(qJ(1));
t238 = t220 * t227;
t223 = sin(qJ(2));
t224 = sin(qJ(1));
t237 = t223 * t224;
t236 = t223 * t227;
t226 = cos(qJ(2));
t235 = t224 * t226;
t234 = t227 * t226;
t233 = qJD(1) * t220;
t222 = sin(qJ(3));
t232 = qJD(2) * t222;
t221 = cos(pkin(6));
t231 = t221 * t234 - t237;
t230 = t221 * t235 + t236;
t229 = t221 * t236 + t235;
t228 = -t221 * t237 + t234;
t219 = t220 * t226 * t232 + (t221 * t222 + t223 * t239) * qJD(3);
t218 = (-t222 * t238 + t229 * t225) * qJD(3) + t231 * t232 + (t228 * t222 - t224 * t239) * qJD(1);
t217 = (t224 * t220 * t222 + t228 * t225) * qJD(3) - t230 * t232 + (-t229 * t222 - t225 * t238) * qJD(1);
t1 = [0, t227 * t233, t231 * qJD(1) + t228 * qJD(2), t217, t217, 0; 0, t224 * t233, t230 * qJD(1) + t229 * qJD(2), t218, t218, 0; 0, 0, t220 * qJD(2) * t223, t219, t219, 0;];
JgD_rot  = t1;
