% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPPR10_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:52
% EndTime: 2019-02-26 22:08:52
% DurationCPUTime: 0.07s
% Computational Cost: add. (9->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
t102 = sin(pkin(6));
t106 = sin(qJ(1));
t115 = t106 * t102;
t105 = sin(qJ(2));
t114 = t106 * t105;
t108 = cos(qJ(2));
t113 = t106 * t108;
t109 = cos(qJ(1));
t112 = t109 * t102;
t111 = t109 * t105;
t110 = t109 * t108;
t107 = cos(qJ(3));
t104 = sin(qJ(3));
t103 = cos(pkin(6));
t1 = [0, t115, t103 * t113 + t111, 0, 0 (-t103 * t114 + t110) * t107 + t104 * t115; 0, -t112, -t103 * t110 + t114, 0, 0 (t103 * t111 + t113) * t107 - t104 * t112; 1, t103, -t102 * t108, 0, 0, t102 * t105 * t107 + t103 * t104;];
Jg_rot  = t1;
