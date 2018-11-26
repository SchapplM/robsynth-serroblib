% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S7RRRRRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% JaD_rot [3x7]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 19:01
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S7RRRRRRR1_jacobiaD_rot_floatb_twist_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_rot_floatb_twist_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_rot_floatb_twist_sym_varpar: qJD has to be [7x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S7RRRRRRR1_jacobiaD_rot_floatb_twist_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiaD_rot_floatb_twist_sym_varpar: pkin has to be [4x1] (double)');
%% Function calls
if link_index == 0
	JaD_rot=S7RRRRRRR1_jacobiaD_rot_0_floatb_twist_sym_varpar(qJ, qJD, pkin);
elseif link_index == 1
	JaD_rot=S7RRRRRRR1_jacobiaD_rot_1_floatb_twist_sym_varpar(qJ, qJD, pkin);
elseif link_index == 2
	JaD_rot=S7RRRRRRR1_jacobiaD_rot_2_floatb_twist_sym_varpar(qJ, qJD, pkin);
elseif link_index == 3
	JaD_rot=S7RRRRRRR1_jacobiaD_rot_3_floatb_twist_sym_varpar(qJ, qJD, pkin);
elseif link_index == 4
	JaD_rot=S7RRRRRRR1_jacobiaD_rot_4_floatb_twist_sym_varpar(qJ, qJD, pkin);
elseif link_index == 5
	JaD_rot=S7RRRRRRR1_jacobiaD_rot_5_floatb_twist_sym_varpar(qJ, qJD, pkin);
elseif link_index == 6
	JaD_rot=S7RRRRRRR1_jacobiaD_rot_6_floatb_twist_sym_varpar(qJ, qJD, pkin);
elseif link_index == 7
	JaD_rot=S7RRRRRRR1_jacobiaD_rot_7_floatb_twist_sym_varpar(qJ, qJD, pkin);
else
	JaD_rot=NaN(6,7);
end