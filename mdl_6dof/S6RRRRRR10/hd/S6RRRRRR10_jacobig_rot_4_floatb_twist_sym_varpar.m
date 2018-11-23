% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR10_jacobig_rot_4_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobig_rot_4_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobig_rot_4_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:18
% EndTime: 2018-11-23 11:27:18
% DurationCPUTime: 0.08s
% Computational Cost: add. (78->27), mult. (87->41), div. (0->0), fcn. (101->19), ass. (0->32)
t238 = sin(pkin(6));
t244 = sin(qJ(1));
t248 = t244 * t238;
t246 = cos(qJ(1));
t247 = t246 * t238;
t245 = cos(qJ(2));
t243 = sin(qJ(2));
t242 = sin(qJ(3));
t241 = cos(pkin(6));
t240 = cos(pkin(7));
t239 = cos(pkin(8));
t237 = sin(pkin(7));
t236 = sin(pkin(8));
t235 = pkin(6) - qJ(2);
t234 = pkin(6) + qJ(2);
t233 = pkin(7) - qJ(3);
t232 = pkin(7) + qJ(3);
t231 = cos(t234);
t230 = sin(t235);
t229 = cos(t235) / 0.2e1;
t228 = sin(t234) / 0.2e1;
t227 = t229 + t231 / 0.2e1;
t226 = cos(t233) / 0.2e1 + cos(t232) / 0.2e1;
t225 = t228 - t230 / 0.2e1;
t224 = t228 + t230 / 0.2e1;
t223 = sin(t232) / 0.2e1 + sin(t233) / 0.2e1;
t222 = -t244 * t227 - t246 * t243;
t221 = t246 * t227 - t244 * t243;
t220 = -t224 * t237 + t241 * t240;
t219 = -t222 * t237 + t240 * t248;
t218 = -t221 * t237 - t240 * t247;
t1 = [0, t248, t219 -(-(-t244 * t225 + t246 * t245) * t242 + t222 * t226 + t223 * t248) * t236 + t219 * t239, 0, 0; 0, -t247, t218 -(-(t246 * t225 + t244 * t245) * t242 + t221 * t226 - t223 * t247) * t236 + t218 * t239, 0, 0; 1, t241, t220 -(-(t229 - t231 / 0.2e1) * t242 + t224 * t226 + t241 * t223) * t236 + t220 * t239, 0, 0;];
Jg_rot  = t1;
