% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPR13_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:37:36
% EndTime: 2019-02-26 22:37:36
% DurationCPUTime: 0.06s
% Computational Cost: add. (38->19), mult. (130->40), div. (0->0), fcn. (134->8), ass. (0->24)
t217 = sin(pkin(6));
t222 = cos(qJ(3));
t236 = t217 * t222;
t224 = cos(qJ(1));
t235 = t217 * t224;
t220 = sin(qJ(2));
t221 = sin(qJ(1));
t234 = t220 * t221;
t233 = t220 * t224;
t223 = cos(qJ(2));
t232 = t221 * t223;
t231 = t224 * t223;
t230 = qJD(1) * t217;
t219 = sin(qJ(3));
t229 = qJD(2) * t219;
t218 = cos(pkin(6));
t228 = t218 * t231 - t234;
t227 = t218 * t232 + t233;
t226 = t218 * t233 + t232;
t225 = -t218 * t234 + t231;
t216 = t217 * t223 * t229 + (t218 * t219 + t220 * t236) * qJD(3);
t215 = (-t219 * t235 + t226 * t222) * qJD(3) + t228 * t229 + (t225 * t219 - t221 * t236) * qJD(1);
t214 = (t221 * t217 * t219 + t225 * t222) * qJD(3) - t227 * t229 + (-t226 * t219 - t222 * t235) * qJD(1);
t1 = [0, t224 * t230, t228 * qJD(1) + t225 * qJD(2), t214, 0, -t214; 0, t221 * t230, t227 * qJD(1) + t226 * qJD(2), t215, 0, -t215; 0, 0, t217 * qJD(2) * t220, t216, 0, -t216;];
JgD_rot  = t1;
