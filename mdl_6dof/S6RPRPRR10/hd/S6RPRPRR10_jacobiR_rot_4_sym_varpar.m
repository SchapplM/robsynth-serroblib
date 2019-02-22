% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:36
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR10_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_jacobiR_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:36:12
% EndTime: 2019-02-22 10:36:12
% DurationCPUTime: 0.03s
% Computational Cost: add. (10->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->11)
t50 = sin(qJ(3));
t51 = sin(qJ(1));
t57 = t51 * t50;
t52 = cos(qJ(3));
t56 = t51 * t52;
t53 = cos(qJ(1));
t55 = t53 * t50;
t54 = t53 * t52;
t49 = cos(pkin(10));
t48 = sin(pkin(10));
t1 = [-t51 * t48 + t49 * t55, 0, t49 * t56, 0, 0, 0; t53 * t48 + t49 * t57, 0, -t49 * t54, 0, 0, 0; 0, 0, -t50 * t49, 0, 0, 0; -t48 * t55 - t51 * t49, 0, -t48 * t56, 0, 0, 0; -t48 * t57 + t53 * t49, 0, t48 * t54, 0, 0, 0; 0, 0, t50 * t48, 0, 0, 0; -t54, 0, t57, 0, 0, 0; -t56, 0, -t55, 0, 0, 0; 0, 0, t52, 0, 0, 0;];
JR_rot  = t1;
