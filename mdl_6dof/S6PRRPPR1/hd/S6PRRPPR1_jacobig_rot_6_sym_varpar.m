% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPPR1_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:13
% EndTime: 2019-02-26 19:58:13
% DurationCPUTime: 0.07s
% Computational Cost: add. (15->10), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->14)
t101 = sin(pkin(10));
t102 = sin(pkin(6));
t110 = t101 * t102;
t103 = cos(pkin(10));
t109 = t103 * t102;
t104 = cos(pkin(6));
t105 = sin(qJ(2));
t108 = t104 * t105;
t106 = cos(qJ(2));
t107 = t104 * t106;
t100 = qJ(3) + pkin(11);
t99 = cos(t100);
t98 = sin(t100);
t1 = [0, t110, t101 * t107 + t103 * t105, 0, 0 (-t101 * t108 + t103 * t106) * t98 - t99 * t110; 0, -t109, t101 * t105 - t103 * t107, 0, 0 (t101 * t106 + t103 * t108) * t98 + t99 * t109; 0, t104, -t102 * t106, 0, 0, t102 * t105 * t98 - t104 * t99;];
Jg_rot  = t1;
