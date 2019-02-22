% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:28
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPPRR3_jacobiR_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobiR_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobiR_rot_3_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:28:44
% EndTime: 2019-02-22 09:28:44
% DurationCPUTime: 0.02s
% Computational Cost: add. (4->4), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
t49 = cos(pkin(6));
t50 = sin(qJ(2));
t53 = t49 * t50;
t51 = cos(qJ(2));
t52 = t49 * t51;
t48 = cos(pkin(10));
t47 = sin(pkin(6));
t46 = sin(pkin(10));
t1 = [0, -t46 * t52 - t48 * t50, 0, 0, 0, 0; 0, -t46 * t50 + t48 * t52, 0, 0, 0, 0; 0, t47 * t51, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, -t46 * t53 + t48 * t51, 0, 0, 0, 0; 0, t46 * t51 + t48 * t53, 0, 0, 0, 0; 0, t47 * t50, 0, 0, 0, 0;];
JR_rot  = t1;
